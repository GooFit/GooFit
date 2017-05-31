/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

TODO:
- Test lineshapes, only done for BW_DP and BW_MINT so far
- Check and implement more SF
- Currently no check if the event is even allowed in phasespace is done. This should preferably be done outside of this
class.

- Some things could be implemented differently maybe, performance should be compared for both cases.
  -For example the way Spinfactors are stored in the same array as the Lineshape values.
   Is this really worth the memory we lose by using a complex to store the SF?
*/

#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include "goofit/Error.h"
#include "goofit/PDFs/physics/DP4Pdf.h"
#include "goofit/PDFs/physics/EvalVar.h"

namespace GooFit {

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!
__device__ thrust::complex<fptype> *cResSF[10];
__device__ thrust::complex<fptype> *Amps_DP[10];
/*
Constant memory array to hold specific info for amplitude calculation.
First entries are the starting points in array, necessary, because number of Lineshapes(LS) or Spinfactors(SF) can vary
|start of each Amplitude| #Linshapes | #Spinfactors | LS-indices | SF-indices|
| 1 entry per Amplitude | 1 per Amp  | 1 per Amp    | #LS in Amp| #SF in Amp|
*/
__constant__ unsigned int AmpIndices[500];

// This function gets called by the GooFit framework to get the value of the PDF.
__device__ fptype device_DP(fptype *evt, fptype *p, unsigned int *indices) {
    // printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real,
    // totalAmp.imag);

    auto evtNum = static_cast<int>(floor(0.5 + evt[indices[7 + indices[0]]]));
    // printf("%i\n",evtNum );
    thrust::complex<fptype> totalAmp(0, 0);
    unsigned int cacheToUse = indices[2];
    unsigned int numAmps    = indices[5];

    for(int i = 0; i < numAmps; ++i) {
        thrust::complex<fptype> amp{p[indices[6 + 2 * i]], p[indices[7 + 2 * i]]};

        thrust::complex<fptype> matrixelement((Amps_DP[cacheToUse][evtNum * numAmps + i]).real(),
                                              (Amps_DP[cacheToUse][evtNum * numAmps + i]).imag());

        totalAmp += matrixelement * amp;
    }

    fptype ret         = thrust::norm(totalAmp);
    int effFunctionIdx = 6 + 2 * indices[3] + 2 * indices[4] + 2 * indices[5];
    fptype eff         = callFunction(evt, indices[effFunctionIdx], indices[effFunctionIdx + 1]);
    ret *= eff;

    // printf("result %.7g\n", ret);
    return ret;
}

__device__ device_function_ptr ptr_to_DP = device_DP;

__host__ DPPdf::DPPdf(std::string n,
                      std::vector<Variable *> observables,
                      DecayInfo_DP *decay,
                      GooPdf *efficiency,
                      unsigned int MCeventsNorm)
    : GooPdf(nullptr, n)
    , decayInfo(decay)
    , _observables(observables)
    , totalEventSize(observables.size()) // number of observables plus eventnumber
{
    for(auto &observable : observables) {
        registerObservable(observable);
    }

    // registerObservable(eventNumber);

    std::vector<fptype> decayConstants;
    decayConstants.push_back(decayInfo->meson_radius);

    for(double &particle_masse : decayInfo->particle_masses) {
        decayConstants.push_back(particle_masse);
    }

    std::vector<unsigned int> pindices;
    pindices.push_back(registerConstants(decayConstants.size()));
    MEMCPY_TO_SYMBOL(functorConstants,
                     &decayConstants[0],
                     decayConstants.size() * sizeof(fptype),
                     cIndex * sizeof(fptype),
                     cudaMemcpyHostToDevice);
    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    pindices.push_back(cacheToUse);
    pindices.push_back(0); //#LS
    pindices.push_back(0); //#SF
    pindices.push_back(0); //#AMP

    // This is the start of reading in the amplitudes and adding the lineshapes and Spinfactors to this PDF
    // This is done in this way so we don't have multiple copies of one lineshape in one pdf.
    std::vector<unsigned int> ampidx;
    std::vector<unsigned int> nPermVec;
    std::vector<unsigned int> ampidxstart;

    for(auto &amplitude : decayInfo->amplitudes) {
        AmpMap[amplitude->_uniqueDecayStr] = std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));

        auto LSvec = amplitude->_LS;

        for(auto &LSIT : LSvec) {
            auto found = std::find_if(components.begin(), components.end(), [&LSIT](const PdfBase *L) {
                return (*LSIT) == *(dynamic_cast<const Lineshape *>(L));
            });

            if(found != components.end()) {
                AmpMap[amplitude->_uniqueDecayStr].first.push_back(std::distance(components.begin(), found));
            } else {
                components.push_back(LSIT);
                AmpMap[amplitude->_uniqueDecayStr].first.push_back(components.size() - 1);
            }
        }

        auto SFvec = amplitude->_SF;

        for(auto &SFIT : SFvec) {
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });

            if(found != SpinFactors.end()) {
                AmpMap[amplitude->_uniqueDecayStr].second.push_back(std::distance(SpinFactors.begin(), found));
            } else {
                SpinFactors.push_back(SFIT);
                AmpMap[amplitude->_uniqueDecayStr].second.push_back(SpinFactors.size() - 1);
            }
        }

        nPermVec.push_back(amplitude->_nPerm);
        pindices.push_back(registerParameter(amplitude->_ar));
        pindices.push_back(registerParameter(amplitude->_ai));

        ampidxstart.push_back(ampidx.size());
        std::vector<unsigned int> ls = AmpMap[amplitude->_uniqueDecayStr].first;
        std::vector<unsigned int> sf = AmpMap[amplitude->_uniqueDecayStr].second;
        ampidx.push_back(ls.size());
        ampidx.push_back(sf.size());
        ampidx.push_back(amplitude->_nPerm);
        ampidx.insert(ampidx.end(), ls.begin(), ls.end());
        ampidx.insert(ampidx.end(), sf.begin(), sf.end());
    }

    MEMCPY_TO_SYMBOL(
        AmpIndices, &(ampidxstart[0]), ampidxstart.size() * sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(AmpIndices,
                     &(ampidx[0]),
                     ampidx.size() * sizeof(unsigned int),
                     ampidxstart.size() * sizeof(unsigned int),
                     cudaMemcpyHostToDevice);

    pindices[2] = components.size();
    pindices[3] = SpinFactors.size();
    pindices[4] = AmpMap.size();

    for(auto &component : components) {
        reinterpret_cast<Lineshape *>(component)->setConstantIndex(cIndex);
        pindices.push_back(reinterpret_cast<Lineshape *>(component)->getFunctionIndex());
        pindices.push_back(reinterpret_cast<Lineshape *>(component)->getParameterIndex());
    }

    for(auto &SpinFactor : SpinFactors) {
        pindices.push_back(SpinFactor->getFunctionIndex());
        pindices.push_back(SpinFactor->getParameterIndex());
        SpinFactor->setConstantIndex(cIndex);
    }

    pindices.push_back(efficiency->getFunctionIndex());
    pindices.push_back(efficiency->getParameterIndex());
    components.push_back(efficiency);

    GET_FUNCTION_ADDR(ptr_to_DP);
    initialize(pindices);

    Integrator   = new NormIntegrator(parameters);
    redoIntegral = new bool[components.size() - 1];
    cachedMasses = new fptype[components.size() - 1];
    cachedWidths = new fptype[components.size() - 1];

    for(int i = 0; i < components.size() - 1; ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        lscalculators.push_back(new LSCalculator(parameters, i));
    }

    for(int i = 0; i < SpinFactors.size(); ++i) {
        sfcalculators.push_back(new SFCalculator(parameters, i));
    }

    for(int i = 0; i < AmpMap.size(); ++i) {
        AmpCalcs.push_back(new AmpCalc(ampidxstart[i], parameters, nPermVec[i]));
    }

    // fprintf(stderr,"#Amp's %i, #LS %i, #SF %i \n", AmpMap.size(), components.size()-1, SpinFactors.size() );

    std::vector<mcbooster::GReal_t> masses(decayInfo->particle_masses.begin() + 1, decayInfo->particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo->particle_masses[0], masses, MCeventsNorm);
    phsp.Generate(mcbooster::Vector4R(decayInfo->particle_masses[0], 0.0, 0.0, 0.0));
    phsp.Unweight();

    auto nAcc                     = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1                       = phsp.GetDaughters(0);
    auto d2                       = phsp.GetDaughters(1);
    auto d3                       = phsp.GetDaughters(2);
    auto d4                       = phsp.GetDaughters(3);

    // auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    // auto zip_end = zip_begin + d1.size();
    // auto new_end = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    printf("After accept-reject we will keep %.i Events for normalization.\n", static_cast<int>(nAcc));
    d1.shrink_to_fit();
    d2.shrink_to_fit();
    d3.shrink_to_fit();
    d4.shrink_to_fit();

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    norm_M12        = mcbooster::RealVector_d(nAcc);
    norm_M34        = mcbooster::RealVector_d(nAcc);
    norm_CosTheta12 = mcbooster::RealVector_d(nAcc);
    norm_CosTheta34 = mcbooster::RealVector_d(nAcc);
    norm_phi        = mcbooster::RealVector_d(nAcc);

    mcbooster::VariableSet_d VarSet(5);
    VarSet[0] = &norm_M12, VarSet[1] = &norm_M34;
    VarSet[2] = &norm_CosTheta12;
    VarSet[3] = &norm_CosTheta34;
    VarSet[4] = &norm_phi;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet);

    norm_SF  = mcbooster::RealVector_d(nAcc * SpinFactors.size());
    norm_LS  = mcbooster::mc_device_vector<thrust::complex<fptype>>(nAcc * (components.size() - 1));
    MCevents = nAcc;

    addSpecialMask(PdfBase::ForceSeparateNorm);
}

// makes the arrays to chache the lineshape values and spinfactors in CachedResSF and the values of the amplitudes in
// cachedAMPs
// I made the choice to have spinfactors necxt to the values of the lineshape in memory. I waste memory by doing this
// because a spinfactor is saved as complex
// It would be nice to test if this is better than having the spinfactors stored seperately.
__host__ void DPPdf::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum = 6
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} should be 3 or more", totalEventSize);

    if(cachedResSF)
        delete cachedResSF;

    if(cachedAMPs)
        delete cachedAMPs;

    numEntries  = dataSize;
    cachedResSF = new thrust::device_vector<thrust::complex<fptype>>(
        dataSize * (components.size() + SpinFactors.size() - 1)); //   -1 because 1 component is efficiency
    void *dummy = thrust::raw_pointer_cast(cachedResSF->data());
    MEMCPY_TO_SYMBOL(cResSF,
                     &dummy,
                     sizeof(thrust::complex<fptype> *),
                     cacheToUse * sizeof(thrust::complex<fptype> *),
                     cudaMemcpyHostToDevice);

    cachedAMPs   = new thrust::device_vector<thrust::complex<fptype>>(dataSize * (AmpCalcs.size()));
    void *dummy2 = thrust::raw_pointer_cast(cachedAMPs->data());
    MEMCPY_TO_SYMBOL(Amps_DP,
                     &dummy2,
                     sizeof(thrust::complex<fptype> *),
                     cacheToUse * sizeof(thrust::complex<fptype> *),
                     cudaMemcpyHostToDevice);

    setForceIntegrals();
}

// this is where the actual magic happens. This function does all the calculations!
__host__ fptype DPPdf::normalize() const {
    // fprintf(stderr, "start normalize\n");
    recursiveSetNormalisation(1); // Not going to normalize efficiency,
    // so set normalisation factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    // check if MINUIT changed any parameters and if so remember that so we know
    // we need to recalculate that lineshape and every amp, that uses that lineshape
    for(unsigned int i = 0; i < components.size() - 1; ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(components[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    SpinsCalculated    = !forceRedoIntegrals;
    forceRedoIntegrals = false;

    // just some thrust iterators for the calculation.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    // Calculate spinfactors only once for normalisation events and real events
    // strided_range is a template implemented in DalitsPlotHelpers.hh
    // it basically goes through the array by increasing the pointer by a certain amount instead of just one step.
    if(!SpinsCalculated) {
        for(int i = 0; i < SpinFactors.size(); ++i) {
            unsigned int offset = components.size() - 1;
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<thrust::complex<fptype>>::iterator>(
                    cachedResSF->begin() + offset + i, cachedResSF->end(), (components.size() + SpinFactors.size() - 1))
                    .begin(),
                *(sfcalculators[i]));

            if(!generation_no_norm) {
                thrust::transform(
                    thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(),
                                                                 norm_M34.begin(),
                                                                 norm_CosTheta12.begin(),
                                                                 norm_CosTheta34.begin(),
                                                                 norm_phi.begin())),
                    thrust::make_zip_iterator(thrust::make_tuple(
                        norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(), norm_CosTheta34.end(), norm_phi.end())),
                    (norm_SF.begin() + (i * MCevents)),
                    NormSpinCalculator(parameters, i));
            }
        }

        SpinsCalculated = true;
    }

    // this calculates the values of the lineshapes and stores them in the array. It is recalculated every time
    // parameters change.
    for(int i = 0; i < components.size() - 1; ++i) {
        if(redoIntegral[i]) {
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<thrust::complex<fptype>>::iterator>(
                    cachedResSF->begin() + i, cachedResSF->end(), (components.size() + SpinFactors.size() - 1))
                    .begin(),
                *(lscalculators[i]));
        }
    }

    // this is a little messy but it basically checks if the amplitude includes one of the recalculated lineshapes and
    // if so recalculates that amplitude
    auto AmpMapIt = AmpMap.begin();

    for(int i = 0; i < AmpCalcs.size(); ++i) {
        std::vector<unsigned int> redoidx((*AmpMapIt).second.first);
        bool redo = false;

        for(unsigned int j : redoidx) {
            if(!redoIntegral[j])
                continue;

            redo = true;
            break;
        }

        if(redo) {
            thrust::transform(eventIndex,
                              eventIndex + numEntries,
                              strided_range<thrust::device_vector<thrust::complex<fptype>>::iterator>(
                                  cachedAMPs->begin() + i, cachedAMPs->end(), AmpCalcs.size())
                                  .begin(),
                              *(AmpCalcs[i]));
        }
    }

    // lineshape value calculation for the normalisation, also recalculated every time parameter change
    if(!generation_no_norm) {
        for(int i = 0; i < components.size() - 1; ++i) {
            if(!redoIntegral[i])
                continue;

            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(),
                                                             norm_M34.begin(),
                                                             norm_CosTheta12.begin(),
                                                             norm_CosTheta34.begin(),
                                                             norm_phi.begin())),
                thrust::make_zip_iterator(thrust::make_tuple(
                    norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(), norm_CosTheta34.end(), norm_phi.end())),
                (norm_LS.begin() + (i * MCevents)),
                NormLSCalculator(parameters, i));
        }
    }

    thrust::constant_iterator<fptype *> normSFaddress(thrust::raw_pointer_cast(norm_SF.data()));
    thrust::constant_iterator<thrust::complex<fptype> *> normLSaddress(thrust::raw_pointer_cast(norm_LS.data()));
    thrust::constant_iterator<int> NumNormEvents(MCevents);

    // this does the rest of the integration with the cached lineshape and spinfactor values for the normalization
    // events
    fptype ret = 1.0;

    if(!generation_no_norm) {
        fptype sumIntegral = 0;
        sumIntegral += thrust::transform_reduce(
            thrust::make_zip_iterator(thrust::make_tuple(eventIndex, NumNormEvents, normSFaddress, normLSaddress)),
            thrust::make_zip_iterator(
                thrust::make_tuple(eventIndex + MCevents, NumNormEvents, normSFaddress, normLSaddress)),
            *Integrator,
            0.,
            thrust::plus<fptype>());
        // MCevents is the number of normalisation events.
        sumIntegral /= MCevents;
        ret = sumIntegral;
    }

    if(std::isnan(ret))
        GooFit::abort(__FILE__, __LINE__, getName() + " NAN normalization in DPPdf", this);

    host_normalisation[parameters] = 1.0 / ret;
    // printf("end of normalize %f\n", ret);
    return ret;
}

__host__
    std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h>
    DPPdf::GenerateSig(unsigned int numEvents) {
    std::vector<mcbooster::GReal_t> masses(decayInfo->particle_masses.begin() + 1, decayInfo->particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo->particle_masses[0], masses, numEvents, generation_offset);
    phsp.Generate(mcbooster::Vector4R(decayInfo->particle_masses[0], 0.0, 0.0, 0.0));

    auto d1 = phsp.GetDaughters(0);
    auto d2 = phsp.GetDaughters(1);
    auto d3 = phsp.GetDaughters(2);
    auto d4 = phsp.GetDaughters(3);

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    auto SigGen_M12_d        = mcbooster::RealVector_d(numEvents);
    auto SigGen_M34_d        = mcbooster::RealVector_d(numEvents);
    auto SigGen_CosTheta12_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_CosTheta34_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_phi_d        = mcbooster::RealVector_d(numEvents);

    mcbooster::VariableSet_d VarSet_d(5);
    VarSet_d[0] = &SigGen_M12_d, VarSet_d[1] = &SigGen_M34_d;
    VarSet_d[2] = &SigGen_CosTheta12_d;
    VarSet_d[3] = &SigGen_CosTheta34_d;
    VarSet_d[4] = &SigGen_phi_d;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet_d);

    auto h1 = new mcbooster::Particles_h(d1);
    auto h2 = new mcbooster::Particles_h(d2);
    auto h3 = new mcbooster::Particles_h(d3);
    auto h4 = new mcbooster::Particles_h(d4);

    mcbooster::ParticlesSet_h ParSet(4);
    ParSet[0] = h1;
    ParSet[1] = h2;
    ParSet[2] = h3;
    ParSet[3] = h4;

    auto SigGen_M12_h        = new mcbooster::RealVector_h(SigGen_M12_d);
    auto SigGen_M34_h        = new mcbooster::RealVector_h(SigGen_M34_d);
    auto SigGen_CosTheta12_h = new mcbooster::RealVector_h(SigGen_CosTheta12_d);
    auto SigGen_CosTheta34_h = new mcbooster::RealVector_h(SigGen_CosTheta34_d);
    auto SigGen_phi_h        = new mcbooster::RealVector_h(SigGen_phi_d);

    mcbooster::VariableSet_h VarSet(5);
    VarSet[0] = SigGen_M12_h, VarSet[1] = SigGen_M34_h;
    VarSet[2] = SigGen_CosTheta12_h;
    VarSet[3] = SigGen_CosTheta34_h;
    VarSet[4] = SigGen_phi_h;

    auto weights = mcbooster::RealVector_d(phsp.GetWeights());
    phsp.~PhaseSpace();

    auto DS = new mcbooster::RealVector_d(6 * numEvents);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 5; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 6);
        thrust::copy(VarSet_d[i]->begin(), VarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 5, DS->end(), 6);
    thrust::copy(eventNumber, eventNumber + numEvents, sr.begin());

    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(numEvents, 6);

    generation_no_norm = true; // we need no normalization for generation, but we do need to make sure that norm = 1;
    SigGenSetIndices();
    copyParams();
    normalize();
    setForceIntegrals();
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    thrust::device_vector<fptype> results(numEvents);
    thrust::constant_iterator<int> eventSize(6);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);

    MetricTaker evalor(this, getMetricPointer("ptr_to_Prob"));
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEvents, arrayAddress, eventSize)),
                      results.begin(),
                      evalor);
    cudaDeviceSynchronize();
    gooFree(dev_event_array);

    thrust::transform(
        results.begin(), results.end(), weights.begin(), weights.begin(), thrust::multiplies<mcbooster::GReal_t>());
    mcbooster::BoolVector_d flags(numEvents);

    thrust::counting_iterator<mcbooster::GLong_t> first(0);
    thrust::counting_iterator<mcbooster::GLong_t> last = first + numEvents;

    auto max = thrust::max_element(weights.begin(), weights.end());
    thrust::transform(first, last, weights.begin(), flags.begin(), mcbooster::FlagAcceptReject((fptype)*max));

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h   = mcbooster::BoolVector_h(flags);
    cudaDeviceSynchronize();
    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}

SFCalculator::SFCalculator(int pIdx, unsigned int sf_idx)
    : _spinfactor_i(sf_idx)
    , _parameters(pIdx) {}

__device__ thrust::complex<fptype> SFCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    int parameter_i       = 6 + (2 * indices[5]) + (indices[3] * 2)
                      + (_spinfactor_i * 2); // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i + 1];

    fptype m12   = evt[indices[2 + indices[0]]];
    fptype m34   = evt[indices[3 + indices[0]]];
    fptype cos12 = evt[indices[4 + indices[0]]];
    fptype cos34 = evt[indices[5 + indices[0]]];
    fptype phi   = evt[indices[6 + indices[0]]];

    fptype vecs[16];
    get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);
    // printf("%i, %i, %f, %f, %f, %f, %f \n",evtNum, thrust::get<2>(t), m12, m34, cos12, cos34, phi );
    // printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
    // printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
    // printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
    // printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);

    auto func = reinterpret_cast<spin_function_ptr>(device_function_table[functn_i]);
    fptype sf = (*func)(vecs, paramIndices + params_i);
    // printf("SpinFactors %i : %.7g\n",evtNum, sf );
    return thrust::complex<fptype>(sf, 0);
}

NormSpinCalculator::NormSpinCalculator(int pIdx, unsigned int sf_idx)
    : _spinfactor_i(sf_idx)
    , _parameters(pIdx) {}

__device__ fptype NormSpinCalculator::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const {
    unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    unsigned int numLS    = indices[3];
    unsigned int numAmps  = indices[5];
    int parameter_i       = 6 + (2 * numAmps) + (numLS * 2)
                      + (_spinfactor_i * 2); // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i + 1];

    fptype m12   = (thrust::get<0>(t));
    fptype m34   = (thrust::get<1>(t));
    fptype cos12 = (thrust::get<2>(t));
    fptype cos34 = (thrust::get<3>(t));
    fptype phi   = (thrust::get<4>(t));

    fptype vecs[16];
    get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);

    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,0, vecs[0], vecs[1], vecs[2], vecs[3]);
    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,1, vecs[4], vecs[5], vecs[6], vecs[7]);
    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,2, vecs[8], vecs[9], vecs[10], vecs[11]);
    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,3, vecs[12], vecs[13], vecs[14], vecs[15]);
    // // }
    auto func = reinterpret_cast<spin_function_ptr>(device_function_table[functn_i]);
    fptype sf = (*func)(vecs, paramIndices + params_i);

    // printf("NormSF evt:%.5g, %.5g, %.5g, %.5g, %.5g\n", m12, m34, cos12, cos34, phi);
    // printf("NormSF %i, %.5g\n",_spinfactor_i, sf );
    // THREAD_SYNCH
    return sf;
}

LSCalculator::LSCalculator(int pIdx, unsigned int res_idx)
    : _resonance_i(res_idx)
    , _parameters(pIdx) {}

__device__ thrust::complex<fptype> LSCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculates the BW values for a specific resonance.
    thrust::complex<fptype> ret;

    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    int parameter_i
        = 6 + (2 * indices[5]) + (_resonance_i * 2); // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i + 1];
    unsigned int pair     = (paramIndices + params_i)[5];

    fptype m1 = functorConstants[indices[1] + 2];
    fptype m2 = functorConstants[indices[1] + 3];
    fptype m3 = functorConstants[indices[1] + 4];
    fptype m4 = functorConstants[indices[1] + 5];

    fptype m12   = evt[indices[2 + indices[0]]];
    fptype m34   = evt[indices[3 + indices[0]]];
    fptype cos12 = evt[indices[4 + indices[0]]];
    fptype cos34 = evt[indices[5 + indices[0]]];
    fptype phi   = evt[indices[6 + indices[0]]];

    if(pair < 2) {
        fptype mres = pair == 0 ? m12 : m34;
        fptype d1   = pair == 0 ? m1 : m3;
        fptype d2   = pair == 0 ? m2 : m4;
        ret         = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
        // printf("LS %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    } else {
        fptype vecs[16];
        get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret         = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
        // printf("LS_m_calc %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    }

    // if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
    // printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    // printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    // printf("%i mass: %.5g, BW_%i : %f %f\n",evtNum, massstore, _resonance_i, ret.real, ret.imag);

    return ret;
}

NormLSCalculator::NormLSCalculator(int pIdx, unsigned int res_idx)
    : _resonance_i(res_idx)
    , _parameters(pIdx) {}

__device__ thrust::complex<fptype> NormLSCalculator::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const {
    // Calculates the BW values for a specific resonance.
    thrust::complex<fptype> ret;

    unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    unsigned int numAmps  = indices[5];
    int parameter_i
        = 6 + (2 * numAmps) + (_resonance_i * 2); // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i + 1];
    unsigned int pair     = (paramIndices + params_i)[5];

    fptype m1 = functorConstants[indices[1] + 2];
    fptype m2 = functorConstants[indices[1] + 3];
    fptype m3 = functorConstants[indices[1] + 4];
    fptype m4 = functorConstants[indices[1] + 5];

    fptype m12   = (thrust::get<0>(t));
    fptype m34   = (thrust::get<1>(t));
    fptype cos12 = (thrust::get<2>(t));
    fptype cos34 = (thrust::get<3>(t));
    fptype phi   = (thrust::get<4>(t));

    if(pair < 2) {
        fptype mres = pair == 0 ? m12 : m34;
        fptype d1   = pair == 0 ? m1 : m3;
        fptype d2   = pair == 0 ? m2 : m4;
        ret         = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
    } else {
        fptype vecs[16];
        get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret         = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
    }

    // printf("NormLS %f, %f, %f, %f, %f \n",m12, m34, cos12, cos34, phi );
    // printf("%i, %i, %i, %i, %i \n",numLS, numSF, numAmps, offset, evtNum );
    // printf("NLS %i, %f, %f\n",_resonance_i,ret.real, ret.imag);

    // printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    // printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    // THREAD_SYNCH
    return ret;
}

AmpCalc::AmpCalc(unsigned int AmpIdx, unsigned int pIdx, unsigned int nPerm)
    : _nPerm(nPerm)
    , _AmpIdx(AmpIdx)
    , _parameters(pIdx) {}

__device__ thrust::complex<fptype> AmpCalc::operator()(thrust::tuple<int, fptype *, int> t) const {
    unsigned int *indices   = paramIndices + _parameters;
    unsigned int cacheToUse = indices[2];
    unsigned int totalLS    = indices[3];
    unsigned int totalSF    = indices[4];
    unsigned int totalAMP   = indices[5];
    unsigned int offset     = totalLS + totalSF;
    unsigned int numLS      = AmpIndices[totalAMP + _AmpIdx];
    unsigned int numSF      = AmpIndices[totalAMP + _AmpIdx + 1];
    unsigned int evtNum     = thrust::get<0>(t);

    thrust::complex<fptype> returnVal(0, 0);
    unsigned int SF_step = numSF / _nPerm;
    unsigned int LS_step = numLS / _nPerm;

    for(int i = 0; i < _nPerm; ++i) {
        thrust::complex<fptype> ret(1, 0);
        thrust::complex<fptype> tmp(1, 0);

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            tmp = (cResSF[cacheToUse][evtNum * offset + AmpIndices[totalAMP + _AmpIdx + 3 + j]]);
            ret *= tmp;
            // printf("Lineshape = (%.7g, %.7g)\n", tmp.real, tmp.imag);
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            tmp = (cResSF[cacheToUse][evtNum * offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 3 + numLS + j]]
                       .real());
            ret *= tmp;
            // printf("SF = (%.7g, %.7g)\n", tmp.real, tmp.imag);
        }

        // printf("Lineshape Product * SF = (%.7g, %.7g)\n", ret.real, ret.imag);

        returnVal += ret;
    }

    returnVal *= (1 / sqrt(static_cast<fptype>(_nPerm)));
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real, returnVal.imag);
    return returnVal;
}

NormIntegrator::NormIntegrator(unsigned int pIdx)
    : _parameters(pIdx) {}

__device__ fptype NormIntegrator::operator()(thrust::tuple<int, int, fptype *, thrust::complex<fptype> *> t) const {
    unsigned int *indices = paramIndices + _parameters;
    unsigned int totalAMP = indices[5];

    unsigned int evtNum             = thrust::get<0>(t);
    unsigned int MCevents           = thrust::get<1>(t);
    fptype *SFnorm                  = thrust::get<2>(t) + evtNum;
    thrust::complex<fptype> *LSnorm = thrust::get<3>(t) + evtNum;

    thrust::complex<fptype> returnVal(0, 0);

    for(int amp = 0; amp < totalAMP; ++amp) {
        unsigned int ampidx  = AmpIndices[amp];
        unsigned int numLS   = AmpIndices[totalAMP + ampidx];
        unsigned int numSF   = AmpIndices[totalAMP + ampidx + 1];
        unsigned int nPerm   = AmpIndices[totalAMP + ampidx + 2];
        unsigned int SF_step = numSF / nPerm;
        unsigned int LS_step = numLS / nPerm;
        thrust::complex<fptype> ret2(0, 0);
        // printf("%i, %i, %i, %i, %i, %i, %i, %i, %i, %f\n",ampidx, amp, numLS, numSF, nPerm,AmpIndices[totalAMP +
        // ampidx + 3 + 0], AmpIndices[totalAMP + ampidx + 3 + 1], AmpIndices[totalAMP + ampidx + 3 + 2],
        // AmpIndices[totalAMP + ampidx + 3 + 3], (1/sqrt((fptype)(nPerm))) );

        for(int j = 0; j < nPerm; ++j) {
            thrust::complex<fptype> ret(1, 0);

            for(int i = j * LS_step; i < (j + 1) * LS_step; ++i) {
                thrust::complex<fptype> matrixelement(LSnorm[AmpIndices[totalAMP + ampidx + 3 + i] * MCevents]);
                // printf("Norm BW %i, %.5g, %.5g\n",AmpIndices[totalAMP + ampidx + 3 + i] , matrixelement.real,
                // matrixelement.imag);
                ret *= matrixelement;
            }

            for(int i = j * SF_step; i < (j + 1) * SF_step; ++i) {
                fptype matrixelement = (SFnorm[AmpIndices[totalAMP + ampidx + 3 + numLS + i] * MCevents]);
                // printf("Norm SF %i, %.5g\n",AmpIndices[totalAMP + ampidx + 3 + i] , matrixelement);
                ret *= matrixelement;
            }

            ret2 += ret;
        }

        thrust::complex<fptype> amp_C{cudaArray[indices[2 * amp + 6]], cudaArray[indices[2 * amp + 7]]};
        ret2 *= (1 / sqrt(static_cast<fptype>(nPerm)));
        // printf("Result Amplitude %i, %.5g, %.5g\n",amp, ret2.real, ret2.imag);
        returnVal += ret2 * amp_C;
    }

    return thrust::norm(returnVal);
}
} // namespace GooFit
