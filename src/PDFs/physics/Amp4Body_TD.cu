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

#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4Body.h>
#include <goofit/PDFs/physics/Amp4Body_TD.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/detail/Complex.h>
#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/Amplitude.h>
#include <goofit/PDFs/physics/detail/AmpCalc_TD.h>
#include <goofit/PDFs/physics/detail/FourDblTupleAdd.h>
#include <goofit/PDFs/physics/detail/LSCalculator_TD.h>
#include <goofit/PDFs/physics/detail/NormIntegrator_TD.h>
#include <goofit/PDFs/physics/detail/NormLSCalculator_TD.h>
#include <goofit/PDFs/physics/detail/SFCalculator_TD.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/utilities/DebugTools.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_DeviceCached.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_HostCached.h>
#include <goofit/MathUtils.h>

#include <cstdarg>
#include <map>

namespace GooFit {

struct genExp {
    fptype gamma;
    unsigned int offset;
    unsigned int seed;

    __host__ __device__ genExp(unsigned int c, fptype d, unsigned int _seed)
        : gamma(d)
        , offset(c)
        , seed(_seed){};

    __host__ __device__ auto operator()(unsigned int x) const -> fptype {
        thrust::random::default_random_engine rand(seed);
        thrust::uniform_real_distribution<fptype> dist(0, 1);

        rand.discard(x + offset);

        return -log(dist(rand)) / gamma;
    }
};

struct exp_functor {
    size_t tmpparam; // index to access decay time in event array
    size_t tmpoff;   // offset for discard() -> should correspond to the number of already generated events/batch size
    fptype gammamin; // effective Gamma
    fptype wmax;     // maximum value for envelope function
    exp_functor(size_t tmpparam, size_t tmpoff, fptype gammamin, fptype wmax)
        : tmpparam(tmpparam)
        , tmpoff(tmpoff)
        , gammamin(gammamin)
        , wmax(wmax) {}

    __device__ auto operator()(thrust::tuple<unsigned int, fptype, fptype *, unsigned int> t) -> bool {
        int evtNum          = thrust::get<0>(t);
        auto val_pdf_to_gen = thrust::get<1>(t);
        fptype *evt         = thrust::get<2>(t) + (evtNum * thrust::get<3>(t));
        /*
        t0: event number
        t1: values of pdf to generate (full model with decay-time dependence)
        t2: pointer to event array
        t3: event size
        */

        fptype time = evt[tmpparam];

        thrust::random::minstd_rand0 rand(1431655765);
        thrust::uniform_real_distribution<fptype> dist(0, 1);
        rand.discard(tmpoff + evtNum);

        // accept-reject method
        auto val_pdf_envelope = exp(-time * gammamin) * wmax;
        auto rand_uniform     = dist(rand);
        if(val_pdf_envelope < thrust::get<1>(t)) {
            printf("ERROR: Amp4Body_TD::exp_functor: value of envelope function smaller than pdf to generate in "
                   "accept-reject method! envelope: %f pdf: %f decay time: %f expo: %f\n",
                   val_pdf_envelope,
                   val_pdf_to_gen,
                   time,
                   exp(-time * gammamin));
        }
        return rand_uniform * val_pdf_envelope < val_pdf_to_gen;
    }
};

void Amp4Body_TD::printSelectedLineshapes(const std::vector<unsigned int> &lsIndices) const {
    std::cout << "Lineshapes:" << std::endl;
    for(int l = 0; l < lsIndices.size(); l++) {
        int index = lsIndices[l];
        std::cout << "LS #" << index << ": " << *_LineShapes[index] << std::endl;
    }
}

void Amp4Body_TD::printSelectedSFs(const std::vector<unsigned int> &sfIndices) const {
    std::cout << "Spin factors:" << std::endl;
    for(int s = 0; s < sfIndices.size(); s++) {
        int index = sfIndices[s];
        std::cout << "SF #" << index << ": " << *(_SpinFactors[index]) << std::endl;
    }
}

void Amp4Body_TD::printAmpMappings() const {
    std::cout << "\nAmplitudes included in Amp4Body_TD model:" << std::endl << std::endl;

    for(int a = 0; a < _AmpCalcs.size(); a++) {
        AmpCalc_TD *ampCalc = _AmpCalcs[a];
        Amplitude *amp      = dynamic_cast<Amplitude *>(components[a]);
        std::cout << "Amplitude # " << a << " (" << amp->_uniqueDecayStr << "):" << std::endl;

        std::vector<unsigned int> lsIndices = ampCalc->getLineshapeIndices(_NUM_AMPLITUDES);
        printSelectedLineshapes(lsIndices);

        std::vector<unsigned int> sfIndices = ampCalc->getSpinFactorIndices(_NUM_AMPLITUDES);
        printSelectedSFs(sfIndices);

        std::cout << std::endl;
    } // end loop over amps
}

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!
__device__ fpcomplex *Amps_TD[10];
/*
Constant memory array to hold specific info for amplitude calculation.
First entries are the starting points in array, necessary, because number of Lineshapes(LS) or Spinfactors(SF) can vary
|start of each Amplitude| #Linshapes | #Spinfactors | LS-indices | SF-indices|
| 1 entry per Amplitude | 1 per Amp  | 1 per Amp    | #LS in Amp| #SF in Amp|
*/
// __constant__ unsigned int AmpIndices_TD[100];

// This function gets called by the GooFit framework to get the value of the PDF.
__device__ auto device_Amp4Body_TD(fptype *evt, ParameterContainer &pc) -> fptype {
    // printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real,
    // totalAmp.imag);

    int id_evt = pc.getObservable(5);

    auto evtNum = static_cast<int>(floor(0.5 + RO_CACHE(evt[id_evt])));
    // GOOFIT_TRACE("Amp4Body_TD: Number of events: {}", evtNum);

    unsigned int cacheToUse = pc.getConstant(5);
    unsigned int numAmps    = pc.getConstant(8);
    unsigned int totalSF_LS = pc.getConstant(10);

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
    fptype _time         = RO_CACHE(evt[id_time]);
    fptype _sigma        = RO_CACHE(evt[id_sigma]);

    AmpA *= _SqWStoRSrate;
    /*printf("%i read time: %.5g x: %.5g y: %.5g \n",evtNum, _time, _xmixing, _ymixing);*/

    fptype term1    = thrust::norm(AmpA) + thrust::norm(AmpB);
    fptype term2    = thrust::norm(AmpA) - thrust::norm(AmpB);
    fpcomplex term3 = AmpA * thrust::conj(AmpB);
    // printf("%i dev %.7g %.7g %.7g %.7g\n", evtNum, norm2(AmpA), norm2(AmpB), term3.real, term3.imag);

    // increment Index
    pc.incrementIndex();

    // increment over all our lineshapes and spinfactors to get to the resolution function
    for(int i = 0; i < totalSF_LS; i++)
        pc.incrementIndex();

    // int effFunctionIdx = 12 + 2 * indices[3] + 2 * indices[4] + 2 * indices[6];
    // int resfctidx      = indices[11];
    // int resfctpar      = effFunctionIdx + 2;

    // resolution function?
    fptype ret = (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
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

__device__ device_function_ptr ptr_to_Amp4Body_TD = device_Amp4Body_TD;

// Build Amp4Body_TD where the MC events used for normalization are stored on the host side
// and where normalization computations are done in batches.
// This is much slower than the other case (where the events used for normalization are stored on the device side)
// but allows you to use a larger # of events for normalization than would otherwise would be possible.
__host__ Amp4Body_TD::Amp4Body_TD(std::string n,
                                  std::vector<Observable> observables,
                                  DecayInfo4t decay,
                                  MixingTimeResolution *Tres,
                                  GooPdf *efficiency,
                                  Observable *mistag,
                                  const std::vector<long> &normSeeds,
                                  unsigned int numNormEventsToGenPerBatch)
    : Amp4Body_TD(
        n,
        observables,
        decay,
        Tres,
        efficiency,
        mistag,
        NormEvents_4Body_HostCached::buildBatches(normSeeds, numNormEventsToGenPerBatch, decay.particle_masses)) {
    GOOFIT_INFO("Built Amp4Body_TD model where the MC events used for normalization are stored on the host side.");
    GOOFIT_INFO("This may result in much longer computation times!");
    GOOFIT_INFO("Use the alternate Amp4Body_TD constructor for maximum speed!");
}

// Build Amp4Body_TD where the MC events used for normalization are stored on the device side.
// This is much faster than the other case (where the events used for normalization are stored on the host side)
// but places a lower limit on the maximum # of events that can be used for normalization.
__host__ Amp4Body_TD::Amp4Body_TD(std::string n,
                                  std::vector<Observable> observables,
                                  DecayInfo4t decay,
                                  MixingTimeResolution *Tres,
                                  GooPdf *efficiency,
                                  Observable *mistag,
                                  long normSeed,
                                  unsigned int numNormEventsToGen)
    : Amp4Body_TD(n,
                  observables,
                  decay,
                  Tres,
                  efficiency,
                  mistag,
                  NormEvents_4Body_DeviceCached::buildBatches({normSeed}, numNormEventsToGen, decay.particle_masses)) {}

// Does common initialization
__host__ Amp4Body_TD::Amp4Body_TD(std::string n,
                                  std::vector<Observable> observables,
                                  DecayInfo4t decay,
                                  MixingTimeResolution *Tres,
                                  GooPdf *efficiency,
                                  Observable *mistag,
                                  const std::vector<NormEvents_4Body_Base *> &normEvents)
    : Amp4BodyBase("Amp4Body_TD", n)
    , _DECAY_INFO(decay)
    , _resolution(Tres)
    , _totalEventSize(observables.size() + 2) // number of observables plus eventnumber
{
    _normEvents.resize(normEvents.size());
    for(int n = 0; n < normEvents.size(); n++) {
        _normEvents[n] = std::unique_ptr<NormEvents_4Body_Base>(normEvents[n]);
    }

    // should include m12, m34, cos12, cos34, phi, eventnumber, dtime, sigmat. In this order!
    for(auto &observable : observables) {
        registerObservable(observable);
    }

    // constantsList.push_back(_DECAY_INFO.meson_radius);

    for(double const &particle_masse : _DECAY_INFO.particle_masses) {
        registerConstant(particle_masse);
    }

    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    registerParameter(_DECAY_INFO._tau);
    registerParameter(_DECAY_INFO._xmixing);
    registerParameter(_DECAY_INFO._ymixing);
    registerParameter(_DECAY_INFO._SqWStoRSrate);

    if(_resolution->getDeviceFunction() < 0)
        throw GooFit::GeneralError("The resolution device function index {} must be more than 0",
                                   _resolution->getDeviceFunction());

    registerConstant(cacheToUse);

    // This is the start of reading in the amplitudes and adding the lineshapes and Spinfactors to this PDF
    // This is done in this way so we don't have multiple copies of one lineshape in one pdf.
    unsigned int coeff_counter = 0;

    std::vector<Amplitude *> AmpBuffer;
    std::vector<Amplitude *> AmpsA = _DECAY_INFO.amplitudes;
    std::vector<Amplitude *> AmpsB = _DECAY_INFO.amplitudes_B;

    std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int>>> AmpMap;

    std::vector<unsigned int> nPermVec;
    std::vector<unsigned int> amp_idx;
    std::vector<unsigned int> amp_idx_start;

    unsigned int total_lineshapes_spinfactors = 0;

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
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                _LineShapes.begin(), _LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });

            if(found != _LineShapes.end()) {
                AmpMap[i->_uniqueDecayStr].first.push_back(std::distance(_LineShapes.begin(), found));
                // printf("LS %s found at %i\n",LSIT->getName().c_str(),std::distance(_LineShapes.begin(), found));
            } else {
                _LineShapes.push_back(LSIT);
                AmpMap[i->_uniqueDecayStr].first.push_back(_LineShapes.size() - 1);
                // printf("Adding LS %s\n",LSIT->getName().c_str());
            }
        }

        auto SFvec = i->_SF;

        for(auto &SFIT : SFvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                _SpinFactors.begin(), _SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });

            if(found != _SpinFactors.end()) {
                AmpMap[i->_uniqueDecayStr].second.push_back(std::distance(_SpinFactors.begin(), found));
            } else {
                _SpinFactors.push_back(SFIT);
                AmpMap[i->_uniqueDecayStr].second.push_back(_SpinFactors.size() - 1);
            }
        }

        unsigned int flag = 0;
        auto inB = std::find_if(AmpsB.begin(), AmpsB.end(), [AmpsA, &i](const Amplitude *A) { return *i == (*A); });

        if(inB != AmpsB.end()) {
            flag = 2;
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
        // AmpBuffer.push_back(i);

        // std::cout << "Amplitude " << i->_uniqueDecayStr << std::endl;
        // printSelectedLineshapes(ls);
        // printSelectedSFs(sf);
        // std::cout << std::endl;
    }

    for(auto &i : AmpsB) {
        unsigned int flag = 1;
        auto inB
            = std::find_if(AmpBuffer.begin(), AmpBuffer.end(), [AmpsB, &i](const Amplitude *A) { return *i == (*A); });

        if(inB != AmpBuffer.end())
            continue;

        AmpMap[i->_uniqueDecayStr] = std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));

        components.push_back(i);
        // AmpBuffer.push_back(i);

        registerParameter(i->_ar);
        registerParameter(i->_ai);

        auto LSvec = i->_LS;

        for(auto &LSIT : LSvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                _LineShapes.begin(), _LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });

            if(found != _LineShapes.end()) {
                AmpMap[i->_uniqueDecayStr].first.push_back(std::distance(_LineShapes.begin(), found));
                // printf("LS %s found at %i\n",LSIT->getName().c_str(),std::distance(_LineShapes.begin(), found));
            } else {
                _LineShapes.push_back(LSIT);
                AmpMap[i->_uniqueDecayStr].first.push_back(_LineShapes.size() - 1);
                // printf("Adding LS %s\n",LSIT->getName().c_str());
            }
        }

        auto SFvec = i->_SF;

        for(auto &SFIT : SFvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                _SpinFactors.begin(), _SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });

            if(found != _SpinFactors.end()) {
                AmpMap[i->_uniqueDecayStr].second.push_back(std::distance(_SpinFactors.begin(), found));
            } else {
                _SpinFactors.push_back(SFIT);
                AmpMap[i->_uniqueDecayStr].second.push_back(_SpinFactors.size() - 1);
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

        // std::cout << "Amplitude " << i->_uniqueDecayStr << std::endl;
        // printSelectedLineshapes(ls);
        // printSelectedSFs(sf);
        // std::cout << std::endl;
    }

    registerConstant(_LineShapes.size());  // #LS
    registerConstant(_SpinFactors.size()); // #SF
    _NUM_AMPLITUDES = components.size();
    registerConstant(_NUM_AMPLITUDES); // # AMP
    registerConstant(coeff_counter); // Number of coefficients, because its not necessary to be equal to number of Amps.
    registerConstant(total_lineshapes_spinfactors);

    components.push_back(_resolution);
    components.push_back(efficiency);

    if(mistag) {
        registerObservable(*mistag);
        _totalEventSize = 9;
        // TODO: This needs to be registered later!
        // registerConstant(1); // Flags existence of mistag
    }

    // In case the resolution function needs parameters, this registers them.
    // resolution->createParameters(pindices, this);

    registerFunction("ptr_to_Amp4Body_TD", ptr_to_Amp4Body_TD);

    initialize();

    for(int i = 0; i < _LineShapes.size(); i++) {
        _lscalculators.push_back(new LSCalculator_TD());
    }

    MEMCPY_TO_SYMBOL(
        AmpIndices, &(amp_idx_start[0]), amp_idx_start.size() * sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(AmpIndices,
                     &(amp_idx[0]),
                     amp_idx.size() * sizeof(unsigned int),
                     amp_idx_start.size() * sizeof(unsigned int),
                     cudaMemcpyHostToDevice);

    for(int i = 0; i < _SpinFactors.size(); ++i) {
        _sfcalculators.push_back(new SFCalculator_TD());
    }

    for(int i = 0; i < _NUM_AMPLITUDES; ++i) {
        _AmpCalcs.push_back(new AmpCalc_TD(nPermVec[i], amp_idx_start[i]));
    }

    // fprintf(stderr,"#Amp's %i, #LS %i, #SF %i \n", AmpMap.size(), components.size()-1, _SpinFactors.size() );

    GOOFIT_INFO("_Lineshapes size (# unique LS): {}", _LineShapes.size());
    GOOFIT_INFO("_SpinFactors size (# unique SF): {}", _SpinFactors.size());
    GOOFIT_INFO("Components size: {}", components.size());
    GOOFIT_INFO("Num. amplitudes: {}", _NUM_AMPLITUDES);
    GOOFIT_INFO("Num. amp. calcs: {}", _AmpCalcs.size());

    setSeparateNorm();

    printAmpMappings();
}

__host__ void Amp4Body_TD::populateArrays() {
    PdfBase::populateArrays();

    // TODO: We need to expand populateArrays so we handle components correctly!
    efficiencyFunction = host_function_table.size() - 1;
}

// makes the arrays to cache the lineshape values and spinfactors in CachedResSF and the values of the amplitudes in
// _cachedAMPs
// I made the choice to have spinfactors necxt to the values of the lineshape in memory. I waste memory by doing this
// because a spinfactor is saved as complex
// It would be nice to test if this is better than having the spinfactors stored separately.
__host__ void Amp4Body_TD::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum = 6
    _totalEventSize = evtSize;
    if(_totalEventSize < 3)
        throw GooFit::GeneralError("_totalEventSize {} must be 3 or more", _totalEventSize);

    if(_cachedResSF)
        delete _cachedResSF;

    if(_cachedAMPs)
        delete _cachedAMPs;

    numEntries   = dataSize;
    _cachedResSF = new thrust::device_vector<fpcomplex>(dataSize * (_LineShapes.size() + _SpinFactors.size()));
    void *dummy  = thrust::raw_pointer_cast(_cachedResSF->data());
    MEMCPY_TO_SYMBOL(cResSF_TD, &dummy, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    GOOFIT_TRACE("cachedAmps size:{}, {}, {}, {}, {}",
                 dataSize,
                 _AmpCalcs.size(),
                 components.size() - 1,
                 _SpinFactors.size(),
                 _LineShapes.size());

    _cachedAMPs  = new thrust::device_vector<fpcomplex>(dataSize * (_AmpCalcs.size()));
    void *dummy2 = thrust::raw_pointer_cast(_cachedAMPs->data());
    MEMCPY_TO_SYMBOL(Amps_TD, &dummy2, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    setForceIntegrals();
}

__host__ void Amp4Body_TD::computeCachedValues(const std::vector<bool> &lineshapeChanged,
                                               const std::vector<bool> &amplitudeComponentChanged) {
    // just some thrust iterators for the calculation.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(_totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    // Calculate spinfactors only once for normalization events and real events
    // strided_range is a template implemented in DalitsPlotHelpers.hh
    // it basically goes through the array by increasing the pointer by a certain amount instead of just one step.
    if(!_SpinsCalculated) {
        for(int i = 0; i < _SpinFactors.size(); ++i) {
            unsigned int offset = _LineShapes.size();
            unsigned int stride = _LineShapes.size() + _SpinFactors.size();

            GOOFIT_TRACE("SpinFactors - stride: {}", stride);
            _sfcalculators[i]->setDalitzId(getFunctionIndex());
            _sfcalculators[i]->setSpinFactorId(_SpinFactors[i]->getFunctionIndex());

            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    _cachedResSF->begin() + offset + i, _cachedResSF->end(), stride)
                    .begin(),
                *(_sfcalculators[i]));

            // DebugTools::printDeviceVecComplexVals(_cachedResSF->begin()+offset+i,
            //       _cachedResSF->end(),
            //       stride,
            //       "Cached res SF (in SF block)",
            //       INT_MAX, 20);
        }
    }

    // this calculates the values of the lineshapes and stores them in the array. It is recalculated every time
    // parameters change.
    for(int i = 0; i < _LineShapes.size(); ++i) {
        if(lineshapeChanged[i]) {
            _lscalculators[i]->setDalitzId(getFunctionIndex());
            _lscalculators[i]->setResonanceId(_LineShapes[i]->getFunctionIndex());

            unsigned int stride = _LineShapes.size() + _SpinFactors.size();

            GOOFIT_TRACE("LineShape[{}] - stride: {}", _LineShapes[i]->getName().c_str(), stride);
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    _cachedResSF->begin() + i, _cachedResSF->end(), stride)
                    .begin(),
                *(_lscalculators[i]));

            // DebugTools::printDeviceVecComplexVals(_cachedResSF->begin()+i,
            //                        _cachedResSF->end(),
            //      stride,
            //      "Cached res SF (in LS block)",
            // INT_MAX, 20);
        }
    }

    // checks if the amplitude includes one of the recalculated lineshapes and if so recalculates that amplitude
    for(int a = 0; a < _NUM_AMPLITUDES; ++a) {
        if(amplitudeComponentChanged[a] == false) {
            continue;
        }

        _AmpCalcs[a]->setDalitzId(getFunctionIndex());

        thrust::transform(eventIndex,
                          eventIndex + numEntries,
                          strided_range<thrust::device_vector<fpcomplex>::iterator>(
                              _cachedAMPs->begin() + a, _cachedAMPs->end(), _AmpCalcs.size())
                              .begin(),
                          *(_AmpCalcs[a]));

        // DebugTools::printDeviceVecComplexVals(_cachedAMPs->begin()+a,
        //      _cachedAMPs->end(),
        //      _AmpCalcs.size(),
        //      "Cached AMPs (in AMP block)",
        // INT_MAX, 20);
    } // end loop over amps
}

__host__ std::vector<bool> Amp4Body_TD::areLineshapesChanged() const {
    // check if MINUIT changed any parameters and if so remember that so we know
    // we need to recalculate that lineshape and every amp, that uses that lineshape
    std::vector<bool> lineshapeChanged(_LineShapes.size(), _forceRedoIntegrals);

    if(!_forceRedoIntegrals) {
        for(int l = 0; l < _LineShapes.size(); l++) {
            if(_LineShapes[l]->parametersChanged()) {
                lineshapeChanged[l] = true;
            }
        }
    }

    return lineshapeChanged;
}

__host__ std::vector<bool> Amp4Body_TD::areAmplitudeComponentsChanged() const {
    // check if MINUIT changed any parameters and if so remember that so we know
    // we need to recalculate that lineshape and every amp, that uses that lineshape
    std::vector<bool> amplitudeComponentChanged(_NUM_AMPLITUDES, _forceRedoIntegrals);

    if(!_forceRedoIntegrals) {
        for(int a = 0; a < _NUM_AMPLITUDES; a++) {
            Amplitude *amp = dynamic_cast<Amplitude *>(components[a]);
            if(amp == NULL) {
                throw GeneralError("Error retrieving amplitude components.");
            }

            if(amp->lineshapeParametersChanged()) {
                amplitudeComponentChanged[a] = true;
            }
        }
    }

    return amplitudeComponentChanged;
}

__host__ int Amp4Body_TD::getNumAccNormEvents() const {
    unsigned int totNumAccNormEvents = 0;
    for(auto const &n : _normEvents) {
        totNumAccNormEvents += n->getNumAccNormEvents();
    }
    return totNumAccNormEvents;
}

__host__ fptype Amp4Body_TD::getSumInitNormEventWeights() const {
    fptype sumWeights = 0;
    for(auto const &n : _normEvents) {
        sumWeights += n->getSumInitWeights();
    }
    return sumWeights;
}

// this is where the actual magic happens. This function does all the calculations!
__host__ auto Amp4Body_TD::normalize() -> fptype {
    if(_cachedResSF == nullptr)
        throw GeneralError("You must call dp.setDataSize(currData.getNumEvents(), N) first!");

    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    host_normalizations.sync(d_normalizations);

    // check if MINUIT changed any parameters and if so remember that so we know
    // we need to recalculate that lineshape and every amp, that uses that lineshape
    std::vector<bool> lineshapeChanged          = areLineshapesChanged();
    std::vector<bool> amplitudeComponentChanged = areAmplitudeComponentsChanged();

    _SpinsCalculated    = !_forceRedoIntegrals;
    _forceRedoIntegrals = false;

    fptype ret = 1.0;

    computeCachedValues(lineshapeChanged, amplitudeComponentChanged);
    if(!_generation_no_norm) {
        bool noCachedNormValuesToCompute
            = _SpinsCalculated
              && std::all_of(lineshapeChanged.cbegin(), lineshapeChanged.cend(), [](bool v) { return !v; });

        fptype tau     = parametersList[0];
        fptype xmixing = parametersList[1];
        fptype ymixing = parametersList[2];

        std::vector<fptype> normResults(_normEvents.size());
        for(int n = 0; n < _normEvents.size(); n++) {
            normResults[n] = _normEvents[n]->computeNorm_TD(noCachedNormValuesToCompute,
                                                            _resolution,
                                                            tau,
                                                            xmixing,
                                                            ymixing,
                                                            getFunctionIndex(),
                                                            _SpinsCalculated,
                                                            lineshapeChanged,
                                                            getSFFunctionIndices(),
                                                            getLSFunctionIndices());
        }
        fptype normResultsSum = MathUtils::doNeumaierSummation(normResults);
        ret                   = normResultsSum / getSumInitNormEventWeights();
    }

    _SpinsCalculated = true;

    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;

    GOOFIT_DEBUG("Norm value: {:.20f}\n", ret);

    return ret;
}

__host__ std::vector<unsigned int> Amp4Body_TD::getSFFunctionIndices() const {
    unsigned int numSF = _SpinFactors.size();

    std::vector<unsigned int> sfFunctionIndices(numSF);

    for(int s = 0; s < numSF; s++) {
        sfFunctionIndices[s] = _SpinFactors[s]->getFunctionIndex();
    }

    return sfFunctionIndices;
}

__host__ std::vector<unsigned int> Amp4Body_TD::getLSFunctionIndices() const {
    unsigned int numLS = _LineShapes.size();

    std::vector<unsigned int> lsFunctionIndices(numLS);

    for(int l = 0; l < numLS; l++) {
        lsFunctionIndices[l] = _LineShapes[l]->getFunctionIndex();
    }

    return lsFunctionIndices;
}

__host__ auto Amp4Body_TD::GenerateSig(unsigned int numEvents, int seed) -> std::
    tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h> {
    initialize();
    copyParams();

    std::vector<mcbooster::GReal_t> masses(_DECAY_INFO.particle_masses.begin() + 1, _DECAY_INFO.particle_masses.end());
    mcbooster::PhaseSpace phsp(_DECAY_INFO.particle_masses[0], masses, numEvents, generation_offset);
    if(seed != 0)
        phsp.SetSeed(seed);
    else
        GOOFIT_INFO("Current generator seed {}, offset {}", phsp.GetSeed(), generation_offset);

    phsp.Generate(mcbooster::Vector4R(_DECAY_INFO.particle_masses[0], 0.0, 0.0, 0.0));

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

    thrust::transform(
        index_sequence_begin, index_sequence_begin + nAcc, dtime_d.begin(), genExp(generation_offset, gammamin, seed));

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

    phsp.FreeResources();
    // QUESTION: why 8 variables per event? Is 8th variable filled?
    // ANSWER: _model_m12, _model_m34,         _model_cos12, _model_cos34,
    //   _model_phi, _model_eventNumber, _model_dtime, _model_sigmat
    // -> 8th parameter is decay-time resolution -> not important when generating with truth resolution
    auto DS = new mcbooster::RealVector_d(8 * nAcc);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 5; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 8);
        thrust::copy(VarSet_d[i]->begin(), VarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 5, DS->end(), 8);
    thrust::copy(eventNumber, eventNumber + nAcc, sr.begin());

    // NOTE/QUESTION: first setting decay time to 0?
    // ANSWER: evaluate envelope function at decay time equal zero
    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr2(DS->begin() + 6, DS->end(), 8);
    thrust::fill_n(sr2.begin(), nAcc, 0);

    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(nAcc, 8);

    _generation_no_norm = true; // we need no normalization for generation, but we do need to make sure that norm = 1;
    SigGenSetIndices();
    normalize();
    setForceIntegrals();
    host_normalizations.sync(d_normalizations);

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

    if(wmax > maxWeight && maxWeight != 0) {
        throw GooFit::GeneralError(
            "WARNING: you just encountered a higher maximum weight than observed in previous iterations.\n"
            "WARNING: Consider recalculating your AccRej flags and acceping based upon these.\n"
            "WARNING: previous weight: {}, new weight: {}\n",
            maxWeight,
            wmax);
    }

    maxWeight = wmax > maxWeight ? wmax : maxWeight;

    // QUESTION: why are we doing this? -> copying decay times, where we evaluate function for accept-reject
    thrust::copy(dtime_d.begin(), dtime_d.end(), sr2.begin());
    dtime_d = mcbooster::RealVector_d();
    thrust::device_vector<fptype> results(nAcc);

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + nAcc, arrayAddress, eventSize)),
                      results.begin(),
                      *logger);

    setFitControl(fc);
    cudaDeviceSynchronize();

    thrust::device_vector<bool> flag2(nAcc); // Should be weight2
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
        flag2.begin(), // Should be weight2
        exp_functor(tmpparam, tmpoff, gammamin, wmax));

    // set weights with weight2

    // Maybe max goes here now

    // Include code from Amp4Body regular here

    gooFree(dev_event_array);

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h   = mcbooster::BoolVector_h(flag2);
    cudaDeviceSynchronize();

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}

} // namespace GooFit
