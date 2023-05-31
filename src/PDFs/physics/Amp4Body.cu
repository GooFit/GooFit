/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficiently tested yet and still under heavy development!

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
#include <goofit/PDFs/physics/Amp4Body.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/Amplitude.h>
#include <goofit/PDFs/physics/detail/AmpCalc.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/PDFs/physics/detail/LSCalculator.h>
#include <goofit/PDFs/physics/detail/NormIntegrator.h>
#include <goofit/PDFs/physics/detail/NormLSCalculator.h>
#include <goofit/PDFs/physics/detail/NormSpinCalculator.h>
#include <goofit/PDFs/physics/detail/SFCalculator.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

#include <cstdarg>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

// This function gets called by the GooFit framework to get the value of the PDF.
__device__ auto device_DP(fptype *evt, ParameterContainer &pc) -> fptype {
    // printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real,
    // totalAmp.imag);

    // TODO: Figure out the offset for the event number observable.
    int id_num  = pc.getObservable(5);
    auto evtNum = static_cast<int>(floor(0.5 + RO_CACHE(evt[id_num])));
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

__host__ Amp4Body::Amp4Body(
    std::string n, std::vector<Observable> observables, DecayInfo4 decay, GooPdf *efficiency, unsigned int MCeventsNorm)
    : Amp4BodyBase("Amp4Body", n)
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
    _NUM_AMPLITUDES = components.size();
    registerConstant(LineShapes.size());            // #LS
    registerConstant(SpinFactors.size());           // #SF
    registerConstant(_NUM_AMPLITUDES);              // #AMP
    registerConstant(total_lineshapes_spinfactors); // total line shapes and spin factors used

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

    for(int i = 0; i < LineShapes.size(); ++i) {
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

    for(int i = 0; i < _NUM_AMPLITUDES; ++i) {
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

    setSeparateNorm();
}

// save our efficiency function.  Resonance's are saved first, then the efficiency function.

__host__ void Amp4Body::populateArrays() {
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
    efficiencyFunction = host_function_table.size() - 1;
}
// makes the arrays to cache the lineshape values and spinfactors in CachedResSF and the values of the amplitudes in
// cachedAMPs
// I made the choice to have spinfactors necxt to the values of the lineshape in memory. I waste memory by doing this
// because a spinfactor is saved as complex
// It would be nice to test if this is better than having the spinfactors stored separately.
__host__ void Amp4Body::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum = 6
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} should be 3 or more", totalEventSize);

    if(cachedResSF)
        delete cachedResSF;

    if(cachedAMPs)
        delete cachedAMPs;

    numEntries = dataSize;

#ifdef GOOFIT_MPI
    int myId, numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    int perTask = numEntries / numProcs;

    int counts[numProcs];
    for(int i = 0; i < numProcs - 1; i++)
        counts[i] = perTask;

    counts[numProcs - 1] = numEntries - perTask * (numProcs - 1);

    setNumPerTask(this, counts[myId]);
#endif

    cachedResSF = new thrust::device_vector<fpcomplex>(
#ifdef GOOFIT_MPI
        m_iEventsPerTask * (LineShapes.size() + SpinFactors.size())); //   -1 because 1 component is efficiency
#else
        dataSize * (LineShapes.size() + SpinFactors.size())); //   -1 because 1 component is efficiency
#endif
    void *dummy = thrust::raw_pointer_cast(cachedResSF->data());
    MEMCPY_TO_SYMBOL(cResSF, &dummy, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

#ifdef GOOFIT_MPI
    cachedAMPs = new thrust::device_vector<fpcomplex>(m_iEventsPerTask * (AmpCalcs.size()));
#else
    cachedAMPs                     = new thrust::device_vector<fpcomplex>(dataSize * (AmpCalcs.size()));
#endif
    void *dummy2 = thrust::raw_pointer_cast(cachedAMPs->data());
    MEMCPY_TO_SYMBOL(Amps_DP, &dummy2, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    setForceIntegrals();
}

// this is where the actual magic happens. This function does all the calculations!
__host__ auto Amp4Body::normalize() -> fptype {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    host_normalizations.sync(d_normalizations);

    // check if MINUIT changed any parameters and if so remember that so we know
    // we need to recalculate that lineshape and every amp, that uses that lineshape
    for(unsigned int i = 0; i < LineShapes.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(LineShapes[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    SpinsCalculated    = !forceRedoIntegrals;
    forceRedoIntegrals = false;

#ifdef GOOFIT_MPI
    unsigned int events_to_process = m_iEventsPerTask;
#else
    unsigned int events_to_process = numEntries;
#endif
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
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + events_to_process, dataArray, eventSize)),
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
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + events_to_process, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedResSF->begin() + i, cachedResSF->end(), stride)
                    .begin(),
                *(lscalculators[i]));
        }
    }
    // this is a little messy but it basically checks if the amplitude includes one of the recalculated lineshapes and
    // if so recalculates that amplitude
    // auto AmpMapIt = AmpMap.begin();
    for(int i = 0; i < _NUM_AMPLITUDES; ++i) {
        if(!redoIntegral[i])
            continue;

        AmpCalcs[i]->setDalitzId(getFunctionIndex());
        // AmpCalcs[i]->setAmplitudeId (amp->getFuncionIndex ());
        thrust::transform(eventIndex,
                          eventIndex + events_to_process,
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
        GooFit::abort(__FILE__, __LINE__, getName() + " NAN normalization in Amp4Body", this);

    host_normalizations.at(normalIdx + 1) = 1.0 / ret;
    cachedNormalization                   = 1.0 / ret;
    return ret;
}

__host__ auto Amp4Body::GenerateSig(unsigned int numEvents, int seed) -> std::
    tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h> {
    // Must configure our functions before any calculations!
    // setupObservables();
    // setIndices();

    initialize();

    std::vector<mcbooster::GReal_t> masses(decayInfo.particle_masses.begin() + 1, decayInfo.particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo.particle_masses[0], masses, numEvents, generation_offset);

    if(seed != 0)
        phsp.SetSeed(seed);
    else
        GOOFIT_INFO("Current generator seed {}, offset {}", phsp.GetSeed(), generation_offset);

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
    host_normalizations.sync(d_normalizations);

    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());

    thrust::device_vector<fptype> results;
    GooPdf::evaluate_with_metric(results);

    ranged_print("Results", results.begin(), results.begin() + 6);

    thrust::transform(
        results.begin(), results.end(), weights.begin(), weights.begin(), thrust::multiplies<mcbooster::GReal_t>());

    mcbooster::BoolVector_d flags(numEvents);
    fillMCFlags(flags, weights, numEvents);

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h   = mcbooster::BoolVector_h(flags);
    cudaDeviceSynchronize();

    setFitControl(fc);

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}

} // namespace GooFit
