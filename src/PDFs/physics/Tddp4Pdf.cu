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
#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/PDFs/physics/EvalVar.h>
#include <goofit/PDFs/physics/Tddp4Pdf.h>
#include <goofit/detail/Complex.h>
#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include <cstdarg>

namespace GooFit {

struct genExp {
    fptype gamma;
    unsigned int offset;

    __host__ __device__ genExp(unsigned int c, fptype d)
        : gamma(d)
        , offset(c){};

    __host__ __device__ fptype operator()(unsigned int x) const {
        thrust::random::default_random_engine rand(1431655765);
        thrust::uniform_real_distribution<fptype> dist(0, 1);

        rand.discard(x + offset);

        return -log(dist(rand)) / gamma;
    }
};

struct exp_functor {
    size_t tmpparam, tmpoff;
    fptype gammamin, wmax;
    exp_functor(size_t tmpparam, size_t tmpoff, fptype gammamin, fptype wmax)
        : tmpparam(tmpparam)
        , tmpoff(tmpoff)
        , gammamin(gammamin)
        , wmax(wmax) {}

    __device__ fptype operator()(thrust::tuple<unsigned int, fptype, fptype *, unsigned int> t) {
        int evtNum  = thrust::get<0>(t);
        fptype *evt = thrust::get<2>(t) + (evtNum * thrust::get<3>(t));
        // unsigned int *indices = paramIndices + tmpparam;

        // while
        //    fptype time = evt[indices[8 + indices[0]]];

        // we are grabbing the time out
        fptype time = evt[tmpparam];

        thrust::random::minstd_rand0 rand(1431655765);
        thrust::uniform_real_distribution<fptype> dist(0, 1);
        rand.discard(tmpoff + evtNum);

        return (dist(rand) * exp(-time * gammamin) * wmax) < thrust::get<1>(t);
    }
};

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!
__device__ fpcomplex *cResSF_TD[10];
__device__ fpcomplex *Amps_TD[10];
/*
Constant memory array to hold specific info for amplitude calculation.
First entries are the starting points in array, necessary, because number of Lineshapes(LS) or Spinfactors(SF) can vary
|start of each Amplitude| #Linshapes | #Spinfactors | LS-indices | SF-indices|
| 1 entry per Amplitude | 1 per Amp  | 1 per Amp    | #LS in Amp| #SF in Amp|
*/
// __constant__ unsigned int AmpIndices_TD[100];

// This function gets called by the GooFit framework to get the value of the PDF.
__device__ fptype device_TDDP4(fptype *evt, ParameterContainer &pc) {
    // printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real,
    // totalAmp.imag);

    int id_evt = pc.getObservable(5);

    auto evtNum = static_cast<int>(floor(0.5 + evt[id_evt]));
    // GOOFIT_TRACE("TDDP4: Number of events: {}", evtNum);

    unsigned int cacheToUse = pc.getConstant(6);
    unsigned int numAmps    = pc.getConstant(9);

    fpcomplex AmpA(0, 0);
    fpcomplex AmpB(0, 0);
    fpcomplex amp_A, amp_B;

    int k = 0;

    for(int i = 0; i < numAmps; ++i) {
        unsigned int start = AmpIndices[i];
        unsigned int flag  = AmpIndices[start + 3 + numAmps];
        fpcomplex temp;

        // printf("flag:%i\n",flag);
        switch(flag) {
        case 0:
            amp_A = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpA += temp * amp_A;
            break;

        case 1:
            amp_B = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpB += temp * amp_B;
            break;

        case 2:
            amp_A = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpA += temp * amp_A;

            ++k;
            amp_B = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpB += temp * amp_B;
            break;
        }
    }

    int id_time  = pc.getObservable(6);
    int id_sigma = pc.getObservable(7);

    fptype _tau          = pc.getParameter(0);
    fptype _xmixing      = pc.getParameter(1);
    fptype _ymixing      = pc.getParameter(2);
    fptype _SqWStoRSrate = pc.getParameter(3);
    fptype _time         = evt[id_time];
    fptype _sigma        = evt[id_sigma];

    AmpA *= _SqWStoRSrate;
    /*printf("%i read time: %.5g x: %.5g y: %.5g \n",evtNum, _time, _xmixing, _ymixing);*/

    fptype term1    = thrust::norm(AmpA) + thrust::norm(AmpB);
    fptype term2    = thrust::norm(AmpA) - thrust::norm(AmpB);
    fpcomplex term3 = AmpA * thrust::conj(AmpB);
    // printf("%i dev %.7g %.7g %.7g %.7g\n", evtNum, norm2(AmpA), norm2(AmpB), term3.real, term3.imag);

    // increment Index
    pc.incrementIndex();

    for(int i = 0; i < numAmps * 4 + numAmps * 2; i++)
        pc.incrementIndex();

    // int effFunctionIdx = 12 + 2 * indices[3] + 2 * indices[4] + 2 * indices[6];
    // int resfctidx      = indices[11];
    // int resfctpar      = effFunctionIdx + 2;

    // resolution function?
    fptype ret = (*(reinterpret_cast<device_resfunction_ptr>(device_function_table[pc.funcIdx])))(
        term1, term2, term3.real(), term3.imag(), _tau, _time, _xmixing, _ymixing, _sigma, pc);

    // increment resolution function
    pc.incrementIndex();

    // efficiency function?
    fptype eff = callFunction(evt, pc);
    /*printf("%i result %.7g, eff %.7g\n",evtNum, ret, eff);*/

    ret *= eff;
    /*printf("in prob: %f\n", ret);*/
    return ret;
}

__device__ device_function_ptr ptr_to_TDDP4 = device_TDDP4;

__host__ TDDP4::TDDP4(std::string n,
                      std::vector<Observable> observables,
                      DecayInfo4t decay,
                      MixingTimeResolution *Tres,
                      GooPdf *efficiency,
                      Observable *mistag,
                      unsigned int MCeventsNorm)
    : GooPdf(n)
    , decayInfo(decay)
    , resolution(Tres)
    , totalEventSize(observables.size() + 2) // number of observables plus eventnumber
{
    // should include m12, m34, cos12, cos34, phi, eventnumber, dtime, sigmat. In this order!
    for(auto &observable : observables) {
        registerObservable(observable);
    }

    constantsList.push_back(decayInfo.meson_radius);

    for(double &particle_masse : decayInfo.particle_masses) {
        registerConstant(particle_masse);
    }

    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    registerParameter(decayInfo._tau);
    registerParameter(decayInfo._xmixing);
    registerParameter(decayInfo._ymixing);
    registerParameter(decayInfo._SqWStoRSrate);

    if(resolution->getDeviceFunction() < 0)
        throw GooFit::GeneralError("The resolution device function index {} must be more than 0",
                                   resolution->getDeviceFunction());

    registerConstant(cacheToUse);
    int lsidx  = registerConstant(0); //#LS
    int sfidx  = registerConstant(0); //#SF
    int ampidx = registerConstant(0); //# AMP
    int coeffidx
        = registerConstant(0); // Number of coefficients, because its not necessary to be equal to number of Amps.

    // This is the start of reading in the amplitudes and adding the lineshapes and Spinfactors to this PDF
    // This is done in this way so we don't have multiple copies of one lineshape in one pdf.
    unsigned int coeff_counter = 0;

    std::vector<Amplitude *> AmpBuffer;
    std::vector<Amplitude *> AmpsA = decayInfo.amplitudes;
    std::vector<Amplitude *> AmpsB = decayInfo.amplitudes_B;

    std::vector<unsigned int> nPermVec;
    std::vector<unsigned int> amp_idx;
    std::vector<unsigned int> amp_idx_start;

    for(auto &i : AmpsA) {
        AmpMap[i->_uniqueDecayStr] = std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));
        AmpBuffer.push_back(i);

        // adding amplitude to components
        components.push_back(i);

        // register the parameters from the amplitudes here
        registerParameter(i->_ar);
        registerParameter(i->_ai);

        auto LSvec = i->_LS;

        for(auto &LSIT : LSvec) {
            auto found = std::find_if(
                LineShapes.begin(), LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });

            if(found != LineShapes.end()) {
                AmpMap[i->_uniqueDecayStr].first.push_back(std::distance(LineShapes.begin(), found));
            } else {
                LineShapes.push_back(LSIT);
                AmpMap[i->_uniqueDecayStr].first.push_back(LineShapes.size() - 1);
            }
        }

        auto SFvec = i->_SF;

        for(auto &SFIT : SFvec) {
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });

            if(found != SpinFactors.end()) {
                AmpMap[i->_uniqueDecayStr].second.push_back(std::distance(SpinFactors.begin(), found));
            } else {
                SpinFactors.push_back(SFIT);
                AmpMap[i->_uniqueDecayStr].second.push_back(SpinFactors.size() - 1);
            }
        }

        unsigned int flag = 0;
        auto inB = std::find_if(AmpsB.begin(), AmpsB.end(), [AmpsA, &i](const Amplitude *A) { return *i == (*A); });

        if(inB != AmpsB.end()) {
            // printf("Found in AmpsB as well: %s\n", (*inB)->_uniqueDecayStr.c_str());
            registerParameter((*inB)->_ar);
            registerParameter((*inB)->_ai);
            ++coeff_counter;
        }

        nPermVec.push_back(i->_nPerm);
        std::vector<unsigned int> ls = AmpMap[i->_uniqueDecayStr].first;
        std::vector<unsigned int> sf = AmpMap[i->_uniqueDecayStr].second;
        amp_idx_start.push_back(amp_idx.size());
        amp_idx.push_back(ls.size());
        amp_idx.push_back(sf.size());
        amp_idx.push_back(i->_nPerm);
        amp_idx.push_back(flag);
        amp_idx.insert(amp_idx.end(), ls.begin(), ls.end());
        amp_idx.insert(amp_idx.end(), sf.begin(), sf.end());
        ++coeff_counter;
        AmpBuffer.push_back(i);
    }

    for(auto &i : AmpsB) {
        unsigned int flag = 1;
        auto inB
            = std::find_if(AmpBuffer.begin(), AmpBuffer.end(), [AmpsB, &i](const Amplitude *A) { return *i == (*A); });

        if(inB != AmpBuffer.end())
            continue;

        AmpMap[i->_uniqueDecayStr] = std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));

        components.push_back(i);
        AmpBuffer.push_back(i);

        registerParameter(i->_ar);
        registerParameter(i->_ai);

        auto LSvec = i->_LS;

        for(auto &LSIT : LSvec) {
            auto found = std::find_if(
                LineShapes.begin(), LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });

            if(found != LineShapes.end()) {
                AmpMap[i->_uniqueDecayStr].first.push_back(std::distance(LineShapes.begin(), found));
            } else {
                LineShapes.push_back(LSIT);
                AmpMap[i->_uniqueDecayStr].first.push_back(LineShapes.size() - 1);
            }
        }

        auto SFvec = i->_SF;

        for(auto &SFIT : SFvec) {
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });

            if(found != SpinFactors.end()) {
                AmpMap[i->_uniqueDecayStr].second.push_back(std::distance(SpinFactors.begin(), found));
            } else {
                SpinFactors.push_back(SFIT);
                AmpMap[i->_uniqueDecayStr].second.push_back(SpinFactors.size() - 1);
            }
        }

        nPermVec.push_back(i->_nPerm);
        ++coeff_counter;
        amp_idx_start.push_back(amp_idx.size());
        std::vector<unsigned int> ls = AmpMap[i->_uniqueDecayStr].first;
        std::vector<unsigned int> sf = AmpMap[i->_uniqueDecayStr].second;
        amp_idx.push_back(ls.size());
        amp_idx.push_back(sf.size());
        amp_idx.push_back(i->_nPerm);
        amp_idx.push_back(flag);
        amp_idx.insert(amp_idx.end(), ls.begin(), ls.end());
        amp_idx.insert(amp_idx.end(), sf.begin(), sf.end());
    }

    constantsList[lsidx]    = LineShapes.size();
    constantsList[sfidx]    = SpinFactors.size();
    constantsList[ampidx]   = components.size();
    constantsList[coeffidx] = coeff_counter;

    components.push_back(resolution);
    components.push_back(efficiency);

    if(mistag) {
        registerObservable(*mistag);
        totalEventSize = 9;
        // TODO: This needs to be registered later!
        // registerConstant(1); // Flags existence of mistag
    }

    // In case the resolution function needs parameters, this registers them.
    // resolution->createParameters(pindices, this);

    initialize();

    Integrator   = new NormIntegrator_TD();
    redoIntegral = new bool[components.size() - 1];
    cachedMasses = new fptype[components.size() - 1];
    cachedWidths = new fptype[components.size() - 1];

    for(auto lineshape : LineShapes) {
        lscalculators.push_back(new LSCalculator_TD());
    }

    //for (int i = 0; i < amp_idx.size(); i++)
    //    printf("%i - %i\n", i, amp_idx[i]);

    MEMCPY_TO_SYMBOL(
        AmpIndices, &(amp_idx_start[0]), amp_idx_start.size() * sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(AmpIndices,
                     &(amp_idx[0]),
                     amp_idx.size() * sizeof(unsigned int),
                     amp_idx_start.size() * sizeof(unsigned int),
                     cudaMemcpyHostToDevice);

    for(int i = 0; i < SpinFactors.size(); ++i) {
        sfcalculators.push_back(new SFCalculator_TD());
    }

    for(int i = 0; i < components.size() - 2; ++i) {
        AmpCalcs.push_back(new AmpCalc_TD(nPermVec[i], amp_idx_start[i]));
    }

    // fprintf(stderr,"#Amp's %i, #LS %i, #SF %i \n", AmpMap.size(), components.size()-1, SpinFactors.size() );

    std::vector<mcbooster::GReal_t> masses(decayInfo.particle_masses.begin() + 1, decayInfo.particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo.particle_masses[0], masses, MCeventsNorm, generation_offset);
    phsp.Generate(mcbooster::Vector4R(decayInfo.particle_masses[0], 0.0, 0.0, 0.0));
    phsp.Unweight();

    auto nAcc                     = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1                       = phsp.GetDaughters(0);
    auto d2                       = phsp.GetDaughters(1);
    auto d3                       = phsp.GetDaughters(2);
    auto d4                       = phsp.GetDaughters(3);

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    auto zip_end   = zip_begin + d1.size();
    auto new_end   = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    d1.erase(thrust::get<0>(new_end.get_iterator_tuple()), d1.end());
    d2.erase(thrust::get<1>(new_end.get_iterator_tuple()), d2.end());
    d3.erase(thrust::get<2>(new_end.get_iterator_tuple()), d3.end());
    d4.erase(thrust::get<3>(new_end.get_iterator_tuple()), d4.end());

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
    norm_LS  = mcbooster::mc_device_vector<fpcomplex>(nAcc * (components.size() - 1));
    MCevents = nAcc;

    addSpecialMask(PdfBase::ForceSeparateNorm);
}

__host__ void TDDP4::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_TDDP4);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_TDDP4");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();

    // go over our amplitudes and actually set index values, update.
    std::vector<unsigned int> amp_idx;
    std::vector<unsigned int> amp_idx_start;

    for(int i = 0; i < components.size() - 2; ++i) {
        auto *amp = dynamic_cast<Amplitude *>(components[i]);

        std::vector<Lineshape *> lineshapes   = amp->getLineShapes();
        std::vector<SpinFactor *> spinfactors = amp->getSpinFactors();

        // printf ("i=%i 0=%i 1=%i 2=%i 3=%i\n", i, amp_idx.size(), lineshapes.size(), spinfactors.size(), amp->_nPerm);
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
                GOOFIT_ERROR("Shouldn't happen, could not find lineshape in array!");
        }

        for(auto &SFIT : spinfactors) {
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });
            if(found != SpinFactors.end()) {
                amp_idx.push_back(std::distance(SpinFactors.begin(), found));
            } else
                GOOFIT_ERROR("Shouldn't happen, could not find spin factor in array!");
        }
    }

    MEMCPY_TO_SYMBOL(AmpIndices,
                     &(amp_idx[0]),
                     amp_idx_start.size() * sizeof(unsigned int),
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
__host__ void TDDP4::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum = 6
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} must be 3 or more", totalEventSize);

    if(cachedResSF)
        delete cachedResSF;

    if(cachedAMPs)
        delete cachedAMPs;

    numEntries  = dataSize;
    cachedResSF = new thrust::device_vector<fpcomplex>(dataSize * (LineShapes.size() + SpinFactors.size()));
    void *dummy = thrust::raw_pointer_cast(cachedResSF->data());
    MEMCPY_TO_SYMBOL(cResSF_TD, &dummy, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    GOOFIT_TRACE("cachedAmps size:{}, {}, {}, {}, {}",
                 dataSize,
                 AmpCalcs.size(),
                 components.size() - 1,
                 SpinFactors.size(),
                 LineShapes.size());

    cachedAMPs   = new thrust::device_vector<fpcomplex>(dataSize * (AmpCalcs.size()));
    void *dummy2 = thrust::raw_pointer_cast(cachedAMPs->data());
    MEMCPY_TO_SYMBOL(Amps_TD, &dummy2, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    setForceIntegrals();
}

// this is where the actual magic happens. This function does all the calculations!
__host__ fptype TDDP4::normalize() const {
    if(cachedResSF == nullptr)
        throw GeneralError("You must call dp.setDataSize(currData.getNumEvents(), N) first!");
    // fprintf(stderr, "start normalize\n");
    recursiveSetNormalisation(1); // Not going to normalize efficiency,
    // so set normalisation factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    MEMCPY_TO_SYMBOL(
        d_normalisations, host_normalisations, totalNormalisations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

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
            unsigned int offset = LineShapes.size();
            unsigned int stride = LineShapes.size() + SpinFactors.size();

            GOOFIT_TRACE("SpinFactors - stride: {}", stride);
            sfcalculators[i]->setDalitzId(getFunctionIndex());
            sfcalculators[i]->setSpinFactorId(SpinFactors[i]->getFunctionIndex());

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

    // fprintf(stderr, "normalize after spins\n");

    // this calculates the values of the lineshapes and stores them in the array. It is recalculated every time
    // parameters change.
    for(int i = 0; i < LineShapes.size(); ++i) {
        if(redoIntegral[i]) {
            lscalculators[i]->setDalitzId(getFunctionIndex());
            lscalculators[i]->setResonanceId(LineShapes[i]->getFunctionIndex());

            unsigned int stride = LineShapes.size() + SpinFactors.size();

            GOOFIT_TRACE("LineShapes - stride: {}", stride);
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedResSF->begin() + i, cachedResSF->end(), stride)
                    .begin(),
                *(lscalculators[i]));
        }
    }

    // fprintf(stderr, "normalize after LS\n");

    // this is a little messy but it basically checks if the amplitude includes one of the recalculated lineshapes and
    // if so recalculates that amplitude
    // auto AmpMapIt = AmpMap.begin();

    for(int i = 0; i < components.size() - 2; ++i) {
        auto *amp = dynamic_cast<Amplitude *>(components[i]);
        // std::vector<unsigned int> redoidx((*AmpMapIt).second.first);

        bool redo = false;
        for(unsigned int j = 0; j < components.size() - 2; j++) {
            if(!redoIntegral[j])
                continue;
            redo = true;
            break;
        }

        if(redo) {
            AmpCalcs[i]->setDalitzId(getFunctionIndex());
            // AmpCalcs[i]->setResonanceId(LineShapes[i]->getFunctionIndex());

            thrust::transform(eventIndex,
                              eventIndex + numEntries,
                              strided_range<thrust::device_vector<fpcomplex>::iterator>(
                                  cachedAMPs->begin() + i, cachedAMPs->end(), AmpCalcs.size())
                                  .begin(),
                              *(AmpCalcs[i]));
        }
    }

    // fprintf(stderr, "normalize after Amps\n");

    // lineshape value calculation for the normalisation, also recalculated every time parameter change
    if(!generation_no_norm) {
        for(int i = 0; i < LineShapes.size(); ++i) {
            if(!redoIntegral[i])
                continue;

            NormLSCalculator_TD ns;
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
    auto ret = 1.0;

    if(!generation_no_norm) {
        thrust::tuple<fptype, fptype, fptype, fptype> dummy(0, 0, 0, 0);
        FourDblTupleAdd MyFourDoubleTupleAdditionFunctor;
        thrust::tuple<fptype, fptype, fptype, fptype> sumIntegral;

        Integrator->setDalitzId(getFunctionIndex());

        sumIntegral = thrust::transform_reduce(
            thrust::make_zip_iterator(thrust::make_tuple(eventIndex, NumNormEvents, normSFaddress, normLSaddress)),
            thrust::make_zip_iterator(
                thrust::make_tuple(eventIndex + MCevents, NumNormEvents, normSFaddress, normLSaddress)),
            *Integrator,
            dummy,
            MyFourDoubleTupleAdditionFunctor);

        GOOFIT_TRACE("sumIntegral={}", sumIntegral);

        // printf("normalize A2/#evts , B2/#evts: %.5g, %.5g\n",thrust::get<0>(sumIntegral)/MCevents,
        // thrust::get<1>(sumIntegral)/MCevents);
        fptype tau     = parametersList[5];
        fptype xmixing = parametersList[8];
        fptype ymixing = parametersList[9];

        ret = resolution->normalisation(thrust::get<0>(sumIntegral),
                                        thrust::get<1>(sumIntegral),
                                        thrust::get<2>(sumIntegral),
                                        thrust::get<3>(sumIntegral),
                                        tau,
                                        xmixing,
                                        ymixing);

        // MCevents is the number of normalisation events.
        ret /= MCevents;
    }

    host_normalisations[parameters] = 1.0 / ret;
    //printf("end of normalize %f\n", ret);
    return ret;
}

__host__
    std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
    TDDP4::GenerateSig(unsigned int numEvents) {
    copyParams();

    std::vector<mcbooster::GReal_t> masses(decayInfo.particle_masses.begin() + 1, decayInfo.particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo.particle_masses[0], masses, numEvents, generation_offset);
    phsp.Generate(mcbooster::Vector4R(decayInfo.particle_masses[0], 0.0, 0.0, 0.0));

    phsp.Unweight();

    auto nAcc                     = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1                       = phsp.GetDaughters(0);
    auto d2                       = phsp.GetDaughters(1);
    auto d3                       = phsp.GetDaughters(2);
    auto d4                       = phsp.GetDaughters(3);

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    auto zip_end   = zip_begin + d1.size();
    auto new_end   = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    flags = mcbooster::BoolVector_d();

    d1.erase(thrust::get<0>(new_end.get_iterator_tuple()), d1.end());
    d2.erase(thrust::get<1>(new_end.get_iterator_tuple()), d2.end());
    d3.erase(thrust::get<2>(new_end.get_iterator_tuple()), d3.end());
    d4.erase(thrust::get<3>(new_end.get_iterator_tuple()), d4.end());

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    auto SigGen_M12_d        = mcbooster::RealVector_d(nAcc);
    auto SigGen_M34_d        = mcbooster::RealVector_d(nAcc);
    auto SigGen_CosTheta12_d = mcbooster::RealVector_d(nAcc);
    auto SigGen_CosTheta34_d = mcbooster::RealVector_d(nAcc);
    auto SigGen_phi_d        = mcbooster::RealVector_d(nAcc);
    auto dtime_d             = mcbooster::RealVector_d(nAcc);

    thrust::counting_iterator<unsigned int> index_sequence_begin(0);

    fptype tau      = parametersList[0].getValue();
    fptype ymixing  = parametersList[2].getValue();
    fptype gammamin = 1.0 / tau - fabs(ymixing) / tau;
    //printf("hostparams: %f, %f\n", tau, ymixing);

    thrust::transform(
        index_sequence_begin, index_sequence_begin + nAcc, dtime_d.begin(), genExp(generation_offset, gammamin));

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
    auto dtime_h             = new mcbooster::RealVector_h(dtime_d);

    mcbooster::VariableSet_h VarSet(6);
    VarSet[0] = SigGen_M12_h;
    VarSet[1] = SigGen_M34_h;
    VarSet[2] = SigGen_CosTheta12_h;
    VarSet[3] = SigGen_CosTheta34_h;
    VarSet[4] = SigGen_phi_h;
    VarSet[5] = dtime_h;

    phsp.~PhaseSpace();

    auto DS = new mcbooster::RealVector_d(8 * nAcc);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 5; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 8);
        thrust::copy(VarSet_d[i]->begin(), VarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 5, DS->end(), 8);
    thrust::copy(eventNumber, eventNumber + nAcc, sr.begin());

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr2(DS->begin() + 6, DS->end(), 8);
    thrust::fill_n(sr2.begin(), nAcc, 0);

    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(nAcc, 8);

    generation_no_norm = true; // we need no normalization for generation, but we do need to make sure that norm = 1;
    SigGenSetIndices();
    normalize();
    setForceIntegrals();
    MEMCPY_TO_SYMBOL(
        d_normalisations, host_normalisations, totalNormalisations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    thrust::device_vector<fptype> weights(nAcc);
    thrust::constant_iterator<int> eventSize(8);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);

    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + nAcc, arrayAddress, eventSize)),
                      weights.begin(),
                      *logger);
    cudaDeviceSynchronize();

    fptype wmax = 1.1 * (fptype)*thrust::max_element(weights.begin(), weights.end());

    if(wmax > maxWeight && maxWeight != 0)
        fprintf(stderr,
                "WARNING: you just encountered a higher maximum weight than observed in previous iterations.\n"
                "WARNING: Consider recalculating your AccRej flags and acceping based upon these.\n"
                "WARNING: previous weight: %.4g, new weight: %.4g\n",
                maxWeight,
                wmax);

    maxWeight = wmax > maxWeight ? wmax : maxWeight;

    thrust::copy(dtime_d.begin(), dtime_d.end(), sr2.begin());

    dtime_d = mcbooster::RealVector_d();
    thrust::device_vector<fptype> results(nAcc);

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + nAcc, arrayAddress, eventSize)),
                      results.begin(),
                      *logger);

    setFitControl(fc);

    cudaDeviceSynchronize();

    thrust::device_vector<fptype> flag2(nAcc);
    thrust::counting_iterator<mcbooster::GLong_t> first(0);
    // thrust::counting_iterator<mcbooster::GLong_t> last = first + nAcc;

    // we do not want to copy the whole class to the GPU so capturing *this is not a great option
    // therefore perpare local copies to capture the variables we need
    unsigned int tmpoff   = generation_offset;
    unsigned int tmpparam = 6;
    wmax                  = maxWeight;

    thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, results.begin(), arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + nAcc, results.end(), arrayAddress, eventSize)),
        flag2.begin(),
        exp_functor(tmpparam, tmpoff, gammamin, wmax));

    gooFree(dev_event_array);

    //printf("Offset: %i und wmax:%.5g\n", generation_offset, wmax);

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h   = mcbooster::BoolVector_h(flag2);
    cudaDeviceSynchronize();

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}

SFCalculator_TD::SFCalculator_TD() = default;

__device__ fpcomplex SFCalculator_TD::operator()(thrust::tuple<int, fptype *, int> t) const {
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    // unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    // int parameter_i       = 12 + (2 * indices[6]) + (indices[3] * 2)
    //                  + (_spinfactor_i * 2); // Find position of this resonance relative to DALITZPLOT start
    // unsigned int functn_i = indices[parameter_i];
    // unsigned int params_i = indices[parameter_i + 1];

    ParameterContainer pc;

    // Increment to TDDP function
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

    fptype M  = pc.getConstant(1);
    fptype m1 = pc.getConstant(2);
    fptype m2 = pc.getConstant(3);
    fptype m3 = pc.getConstant(4);
    fptype m4 = pc.getConstant(5);

    fptype vecs[16];
    get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
    // printf("%i, %i, %f, %f, %f, %f, %f \n",evtNum, thrust::get<2>(t), m12, m34, cos12, cos34, phi );
    // printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
    // printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
    // printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
    // printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);

    // increment to spin factor:
    while(pc.funcIdx < _spinfactor_i)
        pc.incrementIndex();

    auto func = reinterpret_cast<spin_function_ptr>(device_function_table[pc.funcIdx]);
    fptype sf = (*func)(vecs, pc);
    // printf("SpinFactors %i : %.7g\n",_spinfactor_i, sf );
    return {sf, 0.0};
}

NormSpinCalculator_TD::NormSpinCalculator_TD() = default;

__device__ fptype NormSpinCalculator_TD::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const {
    // unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    // int parameter_i       = 12 + (2 * indices[6]) + (indices[3] * 2)
    //                  + (_spinfactor_i * 2); // Find position of this resonance relative to DALITZPLOT start
    // unsigned int functn_i = indices[parameter_i];
    // unsigned int params_i = indices[parameter_i + 1];

    fptype m12   = (thrust::get<0>(t));
    fptype m34   = (thrust::get<1>(t));
    fptype cos12 = (thrust::get<2>(t));
    fptype cos34 = (thrust::get<3>(t));
    fptype phi   = (thrust::get<4>(t));

    ParameterContainer pc;

    fptype M  = pc.getConstant(1);
    fptype m1 = pc.getConstant(2);
    fptype m2 = pc.getConstant(3);
    fptype m3 = pc.getConstant(4);
    fptype m4 = pc.getConstant(5);

    // Increment to TDDP function:
    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

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
    // printf("NormSF %i, %.7g\n",_spinfactor_i, sf );
    // THREAD_SYNCH
    return sf;
}

LSCalculator_TD::LSCalculator_TD() = default;

__device__ fpcomplex LSCalculator_TD::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;

    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    // unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    // int parameter_i
    //    = 12 + (2 * indices[6]) + (_resonance_i * 2); // Find position of this resonance relative to DALITZPLOT start
    // unsigned int functn_i = indices[parameter_i];
    // unsigned int params_i = indices[parameter_i + 1];
    // unsigned int pair     = (paramIndices + params_i)[5];

    ParameterContainer pc;

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

    fptype M  = pc.getConstant(1);
    fptype m1 = pc.getConstant(2);
    fptype m2 = pc.getConstant(3);
    fptype m3 = pc.getConstant(4);
    fptype m4 = pc.getConstant(5);

    while(pc.funcIdx < _resonance_i)
        pc.incrementIndex();

    unsigned int pair = pc.getConstant(2);

    if(pair < 2) {
        fptype mres = pair == 0 ? m12 : m34;
        fptype d1   = pair == 0 ? m1 : m3;
        fptype d2   = pair == 0 ? m2 : m4;
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
        // printf("LS_nt %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    } else {
        fptype vecs[16];
        get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
        // printf("LS %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    }

    // printf("LS %i: %.7g, %.7g, %.7g, %.7g, %.7g \n",evtNum, m12, m34, cos12, cos34, phi );

    // if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
    // printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    // printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    // printf("BW_%i : %f %f\n", _resonance_i, ret.real, ret.imag);

    return ret;
}

NormLSCalculator_TD::NormLSCalculator_TD() = default;

__device__ fpcomplex NormLSCalculator_TD::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;

    // unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    // int parameter_i
    //    = 12 + (2 * indices[6]) + (_resonance_i * 2); // Find position of this resonance relative to DALITZPLOT start
    // unsigned int functn_i = indices[parameter_i];
    // unsigned int params_i = indices[parameter_i + 1];
    // unsigned int pair     = (paramIndices + params_i)[5];

    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    fptype m12   = (thrust::get<0>(t));
    fptype m34   = (thrust::get<1>(t));
    fptype cos12 = (thrust::get<2>(t));
    fptype cos34 = (thrust::get<3>(t));
    fptype phi   = (thrust::get<4>(t));

    fptype M  = pc.getConstant(1);
    fptype m1 = pc.getConstant(2);
    fptype m2 = pc.getConstant(3);
    fptype m3 = pc.getConstant(4);
    fptype m4 = pc.getConstant(5);

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
    // THREAD_SYNCH
    return ret;
}

AmpCalc_TD::AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx)
    : _nPerm(nPerm)
    , _AmpIdx(ampIdx) {}

__device__ fpcomplex AmpCalc_TD::operator()(thrust::tuple<int, fptype *, int> t) const {
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int cacheToUse = pc.getConstant(6);
    unsigned int totalLS    = pc.getConstant(7);
    unsigned int totalSF    = pc.getConstant(8);
    unsigned int totalAMP   = pc.getConstant(9);
    unsigned int offset     = totalLS + totalSF;
    unsigned int numLS      = AmpIndices[totalAMP + _AmpIdx];
    unsigned int numSF      = AmpIndices[totalAMP + _AmpIdx + 1];
    unsigned int evtNum     = thrust::get<0>(t);

    fpcomplex returnVal(0, 0);
    unsigned int SF_step = numSF / _nPerm;
    unsigned int LS_step = numLS / _nPerm;

    for(int i = 0; i < _nPerm; ++i) {
        fpcomplex ret(1, 0);

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 4 + j];
            ret *= (cResSF_TD[cacheToUse][evtNum * offset + idx]);
            // printf("Lineshape %i = (%.7g, %.7g)\n", j, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP +
            // _AmpIdx + 4 + j]]).real, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 +
            // j]]).imag);
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 4 + numLS + j];
            ret *= (cResSF_TD[cacheToUse][evtNum * offset + totalLS + idx].real());
            // printf(" SF = %.7g\n", (cResSF_TD[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 4
            // + numLS + j]].real));
        }

        // printf("Lineshape Product * SF = (%.7g, %.7g)\n", ret.real, ret.imag);

        returnVal += ret;
    }

    returnVal *= (1 / sqrt(static_cast<fptype>(_nPerm)));
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real, returnVal.imag);
    return returnVal;
}

NormIntegrator_TD::NormIntegrator_TD() = default;

__device__ thrust::tuple<fptype, fptype, fptype, fptype> NormIntegrator_TD::
operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const {
    // unsigned int *indices = paramIndices + _parameters;
    // unsigned int totalAMP = indices[5];

    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int totalAMP = pc.getConstant(9);

    unsigned int evtNum   = thrust::get<0>(t);
    unsigned int MCevents = thrust::get<1>(t);
    fptype *SFnorm        = thrust::get<2>(t) + evtNum;
    fpcomplex *LSnorm     = thrust::get<3>(t) + evtNum;

    fpcomplex AmpA(0, 0);
    fpcomplex AmpB(0, 0);
    fpcomplex amp_A, amp_B;

    int k = 0;

    for(int amp = 0; amp < totalAMP; ++amp) {
        unsigned int ampidx  = AmpIndices[amp];
        unsigned int numLS   = AmpIndices[totalAMP + ampidx];
        unsigned int numSF   = AmpIndices[totalAMP + ampidx + 1];
        unsigned int nPerm   = AmpIndices[totalAMP + ampidx + 2];
        unsigned int flag    = AmpIndices[totalAMP + ampidx + 3];
        unsigned int SF_step = numSF / nPerm;
        unsigned int LS_step = numLS / nPerm;
        fpcomplex ret2(0, 0);
        // printf("%i, %i, %i, %i, %i, %i, %i, %i, %i, %f\n",ampidx, amp, numLS, numSF, nPerm,AmpIndices[totalAMP +
        // ampidx + 4 + 0], AmpIndices[totalAMP + ampidx + 4 + 1], AmpIndices[totalAMP + ampidx + 4 + 2],
        // AmpIndices[totalAMP + ampidx + 4 + 3], (1/sqrt((fptype)(nPerm))) );

        for(int j = 0; j < nPerm; ++j) {
            fpcomplex ret(1, 0);

            for(int i = j * LS_step; i < (j + 1) * LS_step; ++i) {
                fpcomplex matrixelement(LSnorm[AmpIndices[totalAMP + ampidx + 3 + i] * MCevents]);
                // printf("Norm BW %i, %.5g, %.5g\n",AmpIndices[totalAMP + ampidx + 4 + i] , matrixelement.real,
                // matrixelement.imag);
                ret *= matrixelement;
            }

            for(int i = j * SF_step; i < (j + 1) * SF_step; ++i) {
                fptype matrixelement = (SFnorm[AmpIndices[totalAMP + ampidx + 3 + numLS + i] * MCevents]);
                // printf("Norm SF %i, %.5g\n",AmpIndices[totalAMP + ampidx + 4 + i] , matrixelement);
                ret *= matrixelement;
            }

            ret2 += ret;
        }

        ret2 *= (1 / sqrt(static_cast<fptype>(nPerm)));
        // printf("Result Amplitude %i, %i, %.5g, %.5g\n",flag, amp, ret2.real, ret2.imag);

        switch(flag) {
        case 0:
            amp_A = fpcomplex(pc.getParameter(12 + 2 * (amp + k)), pc.getParameter(13 + 2 * (amp + k)));
            AmpA += ret2 * amp_A;
            break;

        case 1:
            amp_B = fpcomplex(pc.getParameter(12 + 2 * (amp + k)), pc.getParameter(13 + 2 * (amp + k)));
            AmpB += ret2 * amp_B;
            break;

        case 2:
            amp_A = fpcomplex(pc.getParameter(12 + 2 * (amp + k)), pc.getParameter(13 + 2 * (amp + k)));
            AmpA += ret2 * amp_A;
            ++k;
            amp_B = fpcomplex(pc.getParameter(12 + 2 * (amp + k)), pc.getParameter(13 + 2 * (amp + k)));
            AmpB += ret2 * amp_B;
            break;
        }
    }

    fptype _SqWStoRSrate = pc.getParameter(10);
    AmpA *= _SqWStoRSrate;

    auto AmpAB = AmpA * conj(AmpB);
    return {thrust::norm(AmpA), thrust::norm(AmpB), AmpAB.real(), AmpAB.imag()};
}

} // namespace GooFit
