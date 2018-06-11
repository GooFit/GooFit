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

#include <memory>

#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/PDFs/physics/EvalVar.h>

#include <cstdarg>

namespace GooFit {

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!

__device__ fpcomplex *cResSF[10];
__device__ fpcomplex *Amps_DP[10];
/*
Constant memory array to hold specific info for amplitude calculation.
First entries are the starting points in array, necessary, because number of Lineshapes(LS) or Spinfactors(SF) can vary
|start of each Amplitude| #Linshapes | #Spinfactors | LS-indices | SF-indices|
| 1 entry per Amplitude | 1 per Amp  | 1 per Amp    | #LS in Amp| #SF in Amp|
*/
__constant__ unsigned int AmpIndices[500];

// This function gets called by the GooFit framework to get the value of the PDF.
__device__ fptype device_DP(fptype *evt, ParameterContainer &pc) {
    // printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real,
    // totalAmp.imag);

    // TODO: Figure out the offset for the event number observable.
    int id_num  = pc.getObservable(5);
    auto evtNum = static_cast<int>(floor(0.5 + evt[id_num]));
    // printf("%i\n",evtNum );
    fpcomplex totalAmp(0, 0);
    unsigned int cacheToUse  = pc.getConstant(5);
    unsigned int numAmps     = pc.getConstant(8);
    unsigned int total_LS_SF = pc.getConstant(9);

    for(int i = 0; i < numAmps; ++i) {
        fpcomplex amp{pc.getParameter(2 * i), pc.getParameter(2 * i + 1)};

        fpcomplex matrixelement((Amps_DP[cacheToUse][evtNum * numAmps + i]).real(),
                                (Amps_DP[cacheToUse][evtNum * numAmps + i]).imag());

        totalAmp += matrixelement * amp;
    }

    fptype ret = thrust::norm(totalAmp);

    // TODO: we have increment through all our amps, so no need to find efficiency function

    pc.incrementIndex(1, numAmps * 2, 10, 6, 1);

    // Skip our line shapes and spin factors...
    for(int i = 0; i < total_LS_SF; i++)
        pc.incrementIndex();

    fptype eff = callFunction(evt, pc);
    ret *= eff;

    // printf("result %.7g\n", ret);
    return ret;
}

__device__ device_function_ptr ptr_to_DP = device_DP;

__host__ DPPdf::DPPdf(
    std::string n, std::vector<Observable> observables, DecayInfo4 decay, GooPdf *efficiency, unsigned int MCeventsNorm)
    : GooPdf(n)
    , decayInfo(decay)
    , totalEventSize(observables.size()) // number of observables plus eventnumber
{
    for(auto &observable : observables) {
        registerObservable(observable);
    }

    // registerConstant(decayInfo.meson_radius);

    for(double &particle_masse : decayInfo.particle_masses) {
        registerConstant(particle_masse);
    }

    MEMCPY_TO_SYMBOL(c_meson_radius, &decayInfo.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    static int cacheCount = 0;
    cacheToUse            = cacheCount++;

    registerConstant(cacheToUse);
    int lsidx  = registerConstant(0); //#LS
    int sfidx  = registerConstant(0); //#SF
    int ampidx = registerConstant(0); //#AMP
    int ttlidx = registerConstant(0); // total line shapes and spin factors used

    int total_lineshapes_spinfactors = 0;

    // This is the start of reading in the amplitudes and adding the lineshapes and Spinfactors to this PDF
    // This is done in this way so we don't have multiple copies of one lineshape in one pdf.
    for(auto &amplitude : decayInfo.amplitudes) {
        // printf("uniqueStr:%s\n", amplitude->_uniqueDecayStr.c_str());
        AmpMap[amplitude->_uniqueDecayStr] = std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));

        components.push_back(amplitude);

        // register the parameters from the amplitudes here
        registerParameter(amplitude->_ar);
        registerParameter(amplitude->_ai);

        auto LSvec = amplitude->_LS;

        for(auto &LSIT : LSvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                LineShapes.begin(), LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == *(L); });

            if(found != LineShapes.end()) {
                AmpMap[amplitude->_uniqueDecayStr].first.push_back(std::distance(LineShapes.begin(), found));
            } else {
                // components.push_back (LSIT);
                LineShapes.push_back(LSIT);
                // printf("LS:%s\n", LSIT->getName().c_str());
                AmpMap[amplitude->_uniqueDecayStr].first.push_back(components.size() - 1);
            }
        }

        auto SFvec = amplitude->_SF;

        for(auto &spinfactor : SFvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(SpinFactors.begin(), SpinFactors.end(), [&spinfactor](const SpinFactor *S) {
                return (*spinfactor) == (*S);
            });

            if(found == SpinFactors.end())
                SpinFactors.push_back(spinfactor);
        }
    }

    constantsList[lsidx]  = LineShapes.size();
    constantsList[sfidx]  = SpinFactors.size();
    constantsList[ampidx] = components.size();
    constantsList[ttlidx] = total_lineshapes_spinfactors;

    components.push_back(efficiency);

    registerFunction("ptr_to_DP", ptr_to_DP);

    initialize();

    Integrator   = new NormIntegrator();
    redoIntegral = new bool[components.size() - 1];
    cachedMasses = new fptype[components.size() - 1];
    cachedWidths = new fptype[components.size() - 1];

    std::vector<unsigned int> amp_idx;
    std::vector<unsigned int> amp_idx_start;
    std::vector<unsigned int> nPermVec;

    for(auto &lineshape : LineShapes) {
        lscalculators.push_back(new LSCalculator());
    }

    for(int i = 0; i < components.size() - 1; ++i) {
        auto *amp       = dynamic_cast<Amplitude *>(components[i]);
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;

        std::vector<Lineshape *> lineshapes   = amp->getLineShapes();
        std::vector<SpinFactor *> spinfactors = amp->getSpinFactors();

        // No reason for amp index?
        amp_idx_start.push_back(amp_idx.size());

        amp_idx.push_back(lineshapes.size());
        amp_idx.push_back(spinfactors.size());
        amp_idx.push_back(amp->_nPerm);

        nPermVec.push_back(amp->_nPerm);

        for(auto &LSIT : lineshapes) {
            auto found = std::find_if(
                LineShapes.begin(), LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });
            if(found != LineShapes.end()) {
                amp_idx.push_back(std::distance(LineShapes.begin(), found));
            } else
                printf("Shouldn't happen, could not find lineshape in array!\n");
        }

        for(auto &SFIT : spinfactors) {
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });
            if(found != SpinFactors.end()) {
                amp_idx.push_back(std::distance(SpinFactors.begin(), found));
            } else
                printf("Shouldn't happen, could not find spin factor in array!\n");
        }
    }

    // copy over amp information
    MEMCPY_TO_SYMBOL(
        AmpIndices, &(amp_idx_start[0]), amp_idx_start.size() * sizeof(unsigned int), 0, cudaMemcpyHostToDevice);

    // copy over indexes?
    MEMCPY_TO_SYMBOL(AmpIndices,
                     &(amp_idx[0]),
                     amp_idx.size() * sizeof(unsigned int),
                     amp_idx_start.size() * sizeof(unsigned int),
                     cudaMemcpyHostToDevice);

    for(int i = 0; i < SpinFactors.size(); ++i) {
        sfcalculators.push_back(new SFCalculator());
    }

    for(int i = 0; i < components.size() - 1; ++i) {
        AmpCalcs.push_back(new AmpCalc(nPermVec[i], amp_idx_start[i]));
    }

    // fprintf(stderr,"#Amp's %i, #LS %i, #SF %i \n", AmpMap.size(), components.size()-1, SpinFactors.size() );

    std::vector<mcbooster::GReal_t> masses(decayInfo.particle_masses.begin() + 1, decayInfo.particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo.particle_masses[0], masses, MCeventsNorm);
    phsp.Generate(mcbooster::Vector4R(decayInfo.particle_masses[0], 0.0, 0.0, 0.0));
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
    VarSet[0] = &norm_M12;
    VarSet[1] = &norm_M34;
    VarSet[2] = &norm_CosTheta12;
    VarSet[3] = &norm_CosTheta34;
    VarSet[4] = &norm_phi;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet);

    norm_SF  = mcbooster::RealVector_d(nAcc * SpinFactors.size());
    norm_LS  = mcbooster::mc_device_vector<fpcomplex>(nAcc * (LineShapes.size()));
    MCevents = nAcc;

    addSpecialMask(PdfBase::ForceSeparateNorm);
}

// save our efficiency function.  Resonance's are saved first, then the efficiency function.

__host__ void DPPdf::populateArrays() {
    PdfBase::populateArrays();

    // go over our amplitudes and actually set index values, update.
    std::vector<unsigned int> amp_idx;
    std::vector<unsigned int> amp_idx_start;

    for(int i = 0; i < components.size() - 1; ++i) {
        auto *amp = dynamic_cast<Amplitude *>(components[i]);

        std::vector<Lineshape *> lineshapes   = amp->getLineShapes();
        std::vector<SpinFactor *> spinfactors = amp->getSpinFactors();

        amp_idx_start.push_back(amp_idx.size());

        amp_idx.push_back(lineshapes.size());
        amp_idx.push_back(spinfactors.size());
        amp_idx.push_back(amp->_nPerm);

        for(auto &LSIT : lineshapes) {
            auto found = std::find_if(
                LineShapes.begin(), LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });
            if(found != LineShapes.end()) {
                amp_idx.push_back(std::distance(LineShapes.begin(), found));
            } else
                printf("Shouldn't happen, could not find lineshape in array!\n");
        }

        for(auto &SFIT : spinfactors) {
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });
            if(found != SpinFactors.end()) {
                amp_idx.push_back(std::distance(SpinFactors.begin(), found));
            } else
                printf("Shouldn't happen, could not find spin factor in array!\n");
        }
    }

    MEMCPY_TO_SYMBOL(AmpIndices,
                     &(amp_idx[0]),
                     amp_idx.size() * sizeof(unsigned int),
                     amp_idx_start.size() * sizeof(unsigned int),
                     cudaMemcpyHostToDevice);

    // TODO: We need to expand populateArrays so we handle components correctly!
    efficiencyFunction = num_device_functions - 1;
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
    cachedResSF = new thrust::device_vector<fpcomplex>(
        dataSize * (LineShapes.size() + SpinFactors.size())); //   -1 because 1 component is efficiency
    void *dummy = thrust::raw_pointer_cast(cachedResSF->data());
    MEMCPY_TO_SYMBOL(cResSF, &dummy, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    cachedAMPs   = new thrust::device_vector<fpcomplex>(dataSize * (AmpCalcs.size()));
    void *dummy2 = thrust::raw_pointer_cast(cachedAMPs->data());
    MEMCPY_TO_SYMBOL(Amps_DP, &dummy2, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    setForceIntegrals();
}

// this is where the actual magic happens. This function does all the calculations!
__host__ fptype DPPdf::normalize() {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    // check if MINUIT changed any parameters and if so remember that so we know
    // we need to recalculate that lineshape and every amp, that uses that lineshape
    for(unsigned int i = 0; i < LineShapes.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(LineShapes[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    // SigGenSetIndices();

    SpinsCalculated    = !forceRedoIntegrals;
    forceRedoIntegrals = false;

    // just some thrust iterators for the calculation.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    // Calculate spinfactors only once for normalization events and real events
    // strided_range is a template implemented in DalitsPlotHelpers.hh
    // it basically goes through the array by increasing the pointer by a certain amount instead of just one step.
    if(!SpinsCalculated) {
        for(int i = 0; i < SpinFactors.size(); ++i) {
            unsigned int offset = LineShapes.size();
            sfcalculators[i]->setDalitzId(getFunctionIndex());
            sfcalculators[i]->setSpinFactorId(SpinFactors[i]->getFunctionIndex());
            unsigned int stride = LineShapes.size() + SpinFactors.size();
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedResSF->begin() + offset + i, cachedResSF->end(), stride)
                    .begin(),
                *(sfcalculators[i]));

            if(!generation_no_norm) {
                NormSpinCalculator nsc = NormSpinCalculator();
                nsc.setDalitzId(getFunctionIndex());
                nsc.setSpinFactorId(SpinFactors[i]->getFunctionIndex());
                thrust::transform(
                    thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(),
                                                                 norm_M34.begin(),
                                                                 norm_CosTheta12.begin(),
                                                                 norm_CosTheta34.begin(),
                                                                 norm_phi.begin())),
                    thrust::make_zip_iterator(thrust::make_tuple(
                        norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(), norm_CosTheta34.end(), norm_phi.end())),
                    (norm_SF.begin() + (i * MCevents)),
                    nsc);
            }
        }

        SpinsCalculated = true;
    }

    // this calculates the values of the lineshapes and stores them in the array. It is recalculated every time
    // parameters change.
    for(int i = 0; i < LineShapes.size(); ++i) {
        // auto amp = dynamic_cast<Amplitude*> (components[i]);

        if(redoIntegral[i]) {
            // auto ls = amp->getLineShapes();
            lscalculators[i]->setDalitzId(getFunctionIndex());
            lscalculators[i]->setResonanceId(LineShapes[i]->getFunctionIndex());

            unsigned int stride = LineShapes.size() + SpinFactors.size();

            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedResSF->begin() + i, cachedResSF->end(), stride)
                    .begin(),
                *(lscalculators[i]));
        }
    }

    // this is a little messy but it basically checks if the amplitude includes one of the recalculated lineshapes and
    // if so recalculates that amplitude
    // auto AmpMapIt = AmpMap.begin();

    for(int i = 0; i < LineShapes.size(); ++i) {
        if(!redoIntegral[i])
            continue;

        AmpCalcs[i]->setDalitzId(getFunctionIndex());
        // AmpCalcs[i]->setAmplitudeId (amp->getFuncionIndex ());
        thrust::transform(eventIndex,
                          eventIndex + numEntries,
                          strided_range<thrust::device_vector<fpcomplex>::iterator>(
                              cachedAMPs->begin() + i, cachedAMPs->end(), AmpCalcs.size())
                              .begin(),
                          *(AmpCalcs[i]));
    }

    // lineshape value calculation for the normalization, also recalculated every time parameter change
    if(!generation_no_norm) {
        for(int i = 0; i < LineShapes.size(); ++i) {
            if(!redoIntegral[i])
                continue;

            NormLSCalculator ns;
            ns.setDalitzId(getFunctionIndex());
            ns.setResonanceId(LineShapes[i]->getFunctionIndex());

            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(),
                                                             norm_M34.begin(),
                                                             norm_CosTheta12.begin(),
                                                             norm_CosTheta34.begin(),
                                                             norm_phi.begin())),
                thrust::make_zip_iterator(thrust::make_tuple(
                    norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(), norm_CosTheta34.end(), norm_phi.end())),
                (norm_LS.begin() + (i * MCevents)),
                ns);
        }
    }

    thrust::constant_iterator<fptype *> normSFaddress(thrust::raw_pointer_cast(norm_SF.data()));
    thrust::constant_iterator<fpcomplex *> normLSaddress(thrust::raw_pointer_cast(norm_LS.data()));
    thrust::constant_iterator<int> NumNormEvents(MCevents);

    // this does the rest of the integration with the cached lineshape and spinfactor values for the normalization
    // events
    fptype ret = 1.0;

    if(!generation_no_norm) {
        fptype sumIntegral = 0;
        Integrator->setDalitzId(getFunctionIndex());
        sumIntegral += thrust::transform_reduce(
            thrust::make_zip_iterator(thrust::make_tuple(eventIndex, NumNormEvents, normSFaddress, normLSaddress)),
            thrust::make_zip_iterator(
                thrust::make_tuple(eventIndex + MCevents, NumNormEvents, normSFaddress, normLSaddress)),
            *Integrator,
            0.,
            thrust::plus<fptype>());

        GOOFIT_TRACE("sumIntegral={}", sumIntegral);
        // MCevents is the number of normalization events.
        sumIntegral /= MCevents;
        ret = sumIntegral;
    }

    if(std::isnan(ret))
        GooFit::abort(__FILE__, __LINE__, getName() + " NAN normalization in DPPdf", this);

    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    return ret;
}

__host__
    std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h>
    DPPdf::GenerateSig(unsigned int numEvents) {
    // Must configure our functions before any calculations!
    // setupObservables();
    // setIndices();

    std::vector<mcbooster::GReal_t> masses(decayInfo.particle_masses.begin() + 1, decayInfo.particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo.particle_masses[0], masses, numEvents, generation_offset);
    phsp.Generate(mcbooster::Vector4R(decayInfo.particle_masses[0], 0.0, 0.0, 0.0));

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
    VarSet_d[0] = &SigGen_M12_d;
    VarSet_d[1] = &SigGen_M34_d;
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
    VarSet[0] = SigGen_M12_h;
    VarSet[1] = SigGen_M34_h;
    VarSet[2] = SigGen_CosTheta12_h;
    VarSet[3] = SigGen_CosTheta34_h;
    VarSet[4] = SigGen_phi_h;

    mcbooster::RealVector_d weights(phsp.GetWeights());
    phsp.FreeResources();

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
    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    thrust::device_vector<fptype> results(numEvents);
    thrust::constant_iterator<int> eventSize(6);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);

    // TODO: need to call setIndices (or something) in order to point to ptr_to_Prob, and not ptr_to_Nll
    // MetricTaker evalor(this, getMetricPointer("ptr_to_Prob"));
    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEvents, arrayAddress, eventSize)),
                      results.begin(),
                      *logger);
    cudaDeviceSynchronize();
    gooFree(dev_event_array);

    ranged_print("Results", results.begin(), results.begin() + 6);

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

    setFitControl(fc);

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}

SFCalculator::SFCalculator() = default;

__device__ fpcomplex SFCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    ParameterContainer pc;

    // Increment to DP
    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    int id_m12   = pc.getObservable(0);
    int id_m34   = pc.getObservable(1);
    int id_cos12 = pc.getObservable(2);
    int id_cos34 = pc.getObservable(3);
    int id_phi   = pc.getObservable(4);

    fptype m12   = evt[id_m12];
    fptype m34   = evt[id_m34];
    fptype cos12 = evt[id_cos12];
    fptype cos34 = evt[id_cos34];
    fptype phi   = evt[id_phi];

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);

    fptype vecs[16];
    get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
    // printf("%i, %i, %f, %f, %f, %f, %f \n",evtNum, thrust::get<2>(t), m12, m34, cos12, cos34, phi );
    // printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
    // printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
    // printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
    // printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);

    // loop until our appropriate spin factor
    while(pc.funcIdx < _spinfactor_i)
        pc.incrementIndex();

    auto func = reinterpret_cast<spin_function_ptr>(device_function_table[pc.funcIdx]);
    fptype sf = (*func)(vecs, pc);
    // printf("SpinFactors %i : %.7g\n",evtNum, sf );
    return {sf, 0.0};
}

NormSpinCalculator::NormSpinCalculator() = default;

__device__ fptype NormSpinCalculator::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const {
    fptype m12   = (thrust::get<0>(t));
    fptype m34   = (thrust::get<1>(t));
    fptype cos12 = (thrust::get<2>(t));
    fptype cos34 = (thrust::get<3>(t));
    fptype phi   = (thrust::get<4>(t));

    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);

    fptype vecs[16];
    get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);

    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,0, vecs[0], vecs[1], vecs[2], vecs[3]);
    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,1, vecs[4], vecs[5], vecs[6], vecs[7]);
    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,2, vecs[8], vecs[9], vecs[10], vecs[11]);
    //   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,3, vecs[12], vecs[13], vecs[14], vecs[15]);
    // // }
    while(pc.funcIdx < _spinfactor_i)
        pc.incrementIndex();

    auto func = reinterpret_cast<spin_function_ptr>(device_function_table[pc.funcIdx]);
    fptype sf = (*func)(vecs, pc);

    // printf("NormSF evt:%.5g, %.5g, %.5g, %.5g, %.5g\n", m12, m34, cos12, cos34, phi);
    // printf("NormSF %i, %.5g\n",_spinfactor_i, sf );
    return sf;
}

LSCalculator::LSCalculator() = default;

__device__ fpcomplex LSCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;

    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    ParameterContainer pc;

    // Increment to DP
    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    int id_m12   = pc.getObservable(0);
    int id_m34   = pc.getObservable(1);
    int id_cos12 = pc.getObservable(2);
    int id_cos34 = pc.getObservable(3);
    int id_phi   = pc.getObservable(4);

    fptype m12   = evt[id_m12];
    fptype m34   = evt[id_m34];
    fptype cos12 = evt[id_cos12];
    fptype cos34 = evt[id_cos34];
    fptype phi   = evt[id_phi];

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);

    while(pc.funcIdx < _resonance_i)
        pc.incrementIndex();

    unsigned int pair = pc.getConstant(2);

    if(pair < 2) {
        fptype mres = pair == 0 ? m12 : m34;
        fptype d1   = pair == 0 ? m1 : m3;
        fptype d2   = pair == 0 ? m2 : m4;
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
        // printf("LS %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    } else {
        fptype vecs[16];
        get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
        // printf("LS_m_calc %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    }

    // if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
    // printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    // printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    // printf("%i mass: %.5g, BW_%i : %f %f\n",evtNum, massstore, _resonance_i, ret.real, ret.imag);

    return ret;
}

NormLSCalculator::NormLSCalculator() = default;

__device__ fpcomplex NormLSCalculator::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;

    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    fptype m12   = (thrust::get<0>(t));
    fptype m34   = (thrust::get<1>(t));
    fptype cos12 = (thrust::get<2>(t));
    fptype cos34 = (thrust::get<3>(t));
    fptype phi   = (thrust::get<4>(t));

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);

    // skip to our resonance function
    while(pc.funcIdx < _resonance_i)
        pc.incrementIndex();

    unsigned int pair = pc.getConstant(2);

    if(pair < 2) {
        fptype mres = pair == 0 ? m12 : m34;
        fptype d1   = pair == 0 ? m1 : m3;
        fptype d2   = pair == 0 ? m2 : m4;
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
    } else {
        fptype vecs[16];
        // TODO: What is indices[1]?
        get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
    }

    // printf("NormLS %f, %f, %f, %f, %f \n",m12, m34, cos12, cos34, phi );
    // printf("%i, %i, %i, %i, %i \n",numLS, numSF, numAmps, offset, evtNum );
    // printf("NLS %i, %f, %f\n",_resonance_i,ret.real, ret.imag);

    // printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    // printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    return ret;
}

AmpCalc::AmpCalc(unsigned int nPerm, unsigned int amp)
    : _nPerm(nPerm)
    , _AmpIdx(amp) {}

__device__ fpcomplex AmpCalc::operator()(thrust::tuple<int, fptype *, int> t) const {
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int cacheToUse = pc.getConstant(5);
    unsigned int totalLS    = pc.getConstant(6);
    unsigned int totalSF    = pc.getConstant(7);
    unsigned int totalAMP   = pc.getConstant(8);
    unsigned int offset     = totalLS + totalSF;
    unsigned int numLS      = AmpIndices[totalAMP + _AmpIdx];
    unsigned int numSF      = AmpIndices[totalAMP + _AmpIdx + 1];
    unsigned int evtNum     = thrust::get<0>(t);

    fpcomplex returnVal(0, 0);
    unsigned int SF_step = numSF / _nPerm;
    unsigned int LS_step = numLS / _nPerm;

    for(int i = 0; i < _nPerm; ++i) {
        fpcomplex ret(1, 0);
        fpcomplex tmp(1, 0);

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 3 + j];
            tmp     = (cResSF[cacheToUse][evtNum * offset + idx]);
            ret *= tmp;
            // printf("Lineshape = (%.7g, %.7g)\n", tmp.real, tmp.imag);
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 3 + numLS + j];
            tmp     = (cResSF[cacheToUse][evtNum * offset + totalLS + idx].real());
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

NormIntegrator::NormIntegrator() = default;

__device__ fptype NormIntegrator::operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const {
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int totalAMP = pc.getConstant(8);

    unsigned int evtNum   = thrust::get<0>(t);
    unsigned int MCevents = thrust::get<1>(t);
    fptype *SFnorm        = thrust::get<2>(t) + evtNum;
    fpcomplex *LSnorm     = thrust::get<3>(t) + evtNum;

    fpcomplex returnVal(0, 0);

    for(int amp = 0; amp < totalAMP; ++amp) {
        unsigned int ampidx  = AmpIndices[amp];
        unsigned int numLS   = AmpIndices[totalAMP + ampidx];
        unsigned int numSF   = AmpIndices[totalAMP + ampidx + 1];
        unsigned int nPerm   = AmpIndices[totalAMP + ampidx + 2];
        unsigned int SF_step = numSF / nPerm;
        unsigned int LS_step = numLS / nPerm;
        fpcomplex ret2(0, 0);
        // printf("%i, %i, %i, %i, %i, %i, %i, %i, %i, %f\n",ampidx, amp, numLS, numSF, nPerm,AmpIndices[totalAMP +
        // ampidx + 3 + 0], AmpIndices[totalAMP + ampidx + 3 + 1], AmpIndices[totalAMP + ampidx + 3 + 2],
        // AmpIndices[totalAMP + ampidx + 3 + 3], (1/sqrt((fptype)(nPerm))) );

        for(int j = 0; j < nPerm; ++j) {
            fpcomplex ret(1, 0);

            for(int i = j * LS_step; i < (j + 1) * LS_step; ++i) {
                fpcomplex matrixelement(LSnorm[AmpIndices[totalAMP + ampidx + 3 + i] * MCevents]);
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

        fpcomplex amp_C{pc.getParameter(2 * amp + 0), pc.getParameter(2 * amp + 1)};
        ret2 *= (1 / sqrt(static_cast<fptype>(nPerm)));
        // printf("Result Amplitude %i, %.5g, %.5g\n",amp, ret2.real, ret2.imag);
        returnVal += ret2 * amp_C;
    }

    return thrust::norm(returnVal);
}
} // namespace GooFit
