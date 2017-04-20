/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

TODO:
- Test lineshapes, only done for BW_DP and BW_MINT so far
- Check and implement more SF
- Currently no check if the event is even allowed in phasespace is done. This should preferably be done outside of this class.

- Some things could be implemented differently maybe, performance should be compared for both cases.
  -For example the way Spinfactors are stored in the same array as the Lineshape values.
   Is this really worth the memory we lose by using a complex to store the SF?
*/
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/Generate.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GFunctional.h>
#include "goofit/PDFs/Tddp4Pdf.h"
#include "goofit/PDFs/EvalVar.h"
#include "goofit/PDFs/DP4Pdf.h"


struct genUni {
    fptype low, high;
    unsigned int offset;

    __host__ __device__
    genUni(fptype a, fptype b, unsigned int c) : low(a), high(b), offset(c) {};
    // genUni(fptype a, fptype b) : low(a), high(b) {};

    __host__ __device__
    fptype operator()(unsigned int x) const {

        // unsigned int n = x + 47584732571;
        // n = (n + 0x7ed55d16) + (n << 12);
        // n = (n ^ 0xc761c23c) ^ (n >> 19);
        // n = (n + 0x165667b1) + (n << 5);
        // n = (n + 0xd3a2646c) ^ (n << 9);
        // n = (n + 0xfd7046c5) + (n << 3);
        // n = (n ^ 0xb55a4f09) ^ (n >> 16);
        // thrust::random::default_random_engine rand(1431655765);
        thrust::random::minstd_rand0 rand(1431655765);
        thrust::uniform_real_distribution<fptype> dist(low, high);
        rand.discard(x+offset);
        fptype result = dist(rand);
        // printf("inside gen %u %u %.5g\n",x, offset, result );
        return result;
    }
};



// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!
__device__ devcomplex<fptype>* cResSF_TD[10];
__device__ devcomplex<fptype>* Amps_TD[10];
/*
Constant memory array to hold specific info for amplitude calculation.
First entries are the starting points in array, necessary, because number of Lineshapes(LS) or Spinfactors(SF) can vary
|start of each Amplitude| #Linshapes | #Spinfactors | LS-indices | SF-indices|
| 1 entry per Amplitude | 1 per Amp  | 1 per Amp    | #LS in Amp| #SF in Amp|
*/
// __constant__ unsigned int AmpIndices_TD[100];


// This function gets called by the GooFit framework to get the value of the PDF.
__device__ fptype device_TDDP4(fptype* evt, fptype* p, unsigned int* indices) {
    //printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real, totalAmp.imag);

    int evtNum = (int) FLOOR(0.5 + evt[indices[7 + indices[0]]]);
    // printf("%i\n",evtNum );
    unsigned int cacheToUse    = indices[2];
    unsigned int numAmps       = indices[5];

    devcomplex<fptype> AmpA(0, 0);
    devcomplex<fptype> AmpB(0, 0);
    fptype amp_real_A, amp_imag_A, amp_real_B, amp_imag_B;

    int k = 0;

    for(int i = 0; i < numAmps; ++i) {
        unsigned int start = AmpIndices[i];
        unsigned int flag = AmpIndices[start + 3 + numAmps];
        devcomplex<fptype> temp;

        // printf("flag:%i\n",flag);
        switch(flag) {
        case 0:
            amp_real_A = p[indices[12 + 2*(i+k)]];
            amp_imag_A = p[indices[13 + 2*(i+k)]];
            temp = Amps_TD[cacheToUse][evtNum*numAmps + i];
            AmpA += (temp.multiply(amp_real_A, amp_imag_A));
            // printf("DEV0: %.5g, %.5g, %.5g +i %.5g \n", amp_real_A, amp_imag_A, Amps_TD[cacheToUse][evtNum*numAmps + i].real, Amps_TD[cacheToUse][evtNum*numAmps + i].imag);
            break;

        case 1:
            amp_real_B = p[indices[12 + 2*(i+k)]];
            amp_imag_B = p[indices[13 + 2*(i+k)]];
            temp = Amps_TD[cacheToUse][evtNum*numAmps + i];
            AmpB += (temp.multiply(amp_real_B, amp_imag_B));
            // printf("DEV1: %.5g, %.5g, %.5g +i %.5g \n", amp_real_B, amp_imag_B, Amps_TD[cacheToUse][evtNum*numAmps + i].real, Amps_TD[cacheToUse][evtNum*numAmps + i].imag);
            break;

        case 2:
            amp_real_A = p[indices[12 + 2*(i+k)]];
            amp_imag_A = p[indices[13 + 2*(i+k)]];
            temp = Amps_TD[cacheToUse][evtNum*numAmps + i];
            AmpA += (temp.multiply(amp_real_A, amp_imag_A));
            // printf("DEV2.1: %.5g, %.5g, %.5g +i %.5g \n", amp_real_A, amp_imag_A, Amps_TD[cacheToUse][evtNum*numAmps + i].real, Amps_TD[cacheToUse][evtNum*numAmps + i].imag);
            ++k;
            amp_real_B = p[indices[12 + 2*(i+k)]];
            amp_imag_B = p[indices[13 + 2*(i+k)]];
            temp = Amps_TD[cacheToUse][evtNum*numAmps + i];
            AmpB += (temp.multiply(amp_real_B, amp_imag_B));
            // printf("DEV2.2: %.5g, %.5g, %.5g +i %.5g \n", amp_real_B, amp_imag_B, Amps_TD[cacheToUse][evtNum*numAmps + i].real, Amps_TD[cacheToUse][evtNum*numAmps + i].imag);
            break;
        }
    }

    fptype _tau     = p[indices[7]];
    fptype _xmixing = p[indices[8]];
    fptype _ymixing = p[indices[9]];
    fptype _SqWStoRSrate = p[indices[10]];
    fptype _time    = evt[indices[8 + indices[0]]];
    fptype _sigma   = evt[indices[9 + indices[0]]];

    AmpA *= _SqWStoRSrate;
    // printf("%i read time: %.5g x: %.5g y: %.5g \n",evtNum, _time, _xmixing, _ymixing);

    fptype term1    = norm2(AmpA) + norm2(AmpB);
    fptype term2    = norm2(AmpA) - norm2(AmpB);
    devcomplex<fptype> term3    = AmpA * conj(AmpB);
    // printf("%i dev %.7g %.7g %.7g %.7g\n", evtNum, norm2(AmpA), norm2(AmpB), term3.real, term3.imag);

    int effFunctionIdx = 12 + 2*indices[3] + 2*indices[4] + 2*indices[6];
    int resfctidx = indices[11];
    int resfctpar = effFunctionIdx + 2;


    fptype ret = (*(reinterpret_cast<device_resfunction_ptr>(device_function_table[resfctidx])))(term1, term2, term3.real,
                 term3.imag,
                 _tau, _time, _xmixing, _ymixing, _sigma,
                 p, indices + resfctpar);
    fptype eff = callFunction(evt, indices[effFunctionIdx], indices[effFunctionIdx + 1]);
    // printf("%i result %.7g, eff %.7g\n",evtNum, ret, eff);
    ret *= eff;

    return ret;
}

__device__ device_function_ptr ptr_to_TDDP4 = device_TDDP4;

__host__ TDDP4::TDDP4(std::string n,
                      std::vector<Variable*> observables,
                      DecayInfo_DP* decay,
                      MixingTimeResolution* Tres,
                      GooPdf* efficiency,
                      Variable* mistag,
                      unsigned int MCeventsNorm)
    : GooPdf(0, n)
    , decayInfo(decay)
    , _observables(observables)
    , cachedAMPs(0)
    , cachedResSF(0)
    , forceRedoIntegrals(true)
    , totalEventSize(observables.size()+2) // number of observables plus eventnumber
    , cacheToUse(0)
    , SpinsCalculated(false)
    , resolution(Tres)
    , generation_offset(25031992)
    , genlow(0)
    , genhigh(5)
    , generation_no_norm(false) {
    // should include m12, m34, cos12, cos34, phi, eventnumber, dtime, sigmat. In this order!
    for(std::vector<Variable*>::iterator obsIT = observables.begin(); obsIT != observables.end(); ++obsIT) {
        registerObservable(*obsIT);
    }


    std::vector<fptype>decayConstants;
    decayConstants.push_back(decayInfo->meson_radius);

    for(std::vector<fptype>::iterator pmIT = decayInfo->particle_masses.begin(); pmIT != decayInfo->particle_masses.end();
            ++pmIT) {
        decayConstants.push_back(*pmIT);
    }

    if(mistag) {
        registerObservable(mistag);
        totalEventSize = 9;
        decayConstants.push_back(1); // Flags existence of mistag
    }

    std::vector<unsigned int> pindices;
    pindices.push_back(registerConstants(decayConstants.size()));
    MEMCPY_TO_SYMBOL(functorConstants, &decayConstants[0], decayConstants.size()*sizeof(fptype), cIndex*sizeof(fptype),
                     cudaMemcpyHostToDevice);
    static int cacheCount = 0;
    cacheToUse = cacheCount++;
    pindices.push_back(cacheToUse);
    pindices.push_back(0); //#LS
    pindices.push_back(0); //#SF
    pindices.push_back(0); //#AMP
    pindices.push_back(0); //number of coefficients, because its not necessary to be equal to number of Amps.
    pindices.push_back(registerParameter(decayInfo->_tau));
    pindices.push_back(registerParameter(decayInfo->_xmixing));
    pindices.push_back(registerParameter(decayInfo->_ymixing));
    pindices.push_back(registerParameter(decayInfo->_SqWStoRSrate));
    assert(resolution->getDeviceFunction() >= 0);
    pindices.push_back((unsigned int) resolution->getDeviceFunction());

    // This is the start of reading in the amplitudes and adding the lineshapes and Spinfactors to this PDF
    // This is done in this way so we don't have multiple copies of one lineshape in one pdf.
    std::vector<unsigned int> ampidx;
    std::vector<unsigned int> nPermVec;
    std::vector<unsigned int> ampidxstart;
    unsigned int coeff_counter=0;
    std::vector<Amplitude*> AmpBuffer;

    std::vector<Amplitude*> AmpsA = decayInfo->amplitudes;
    std::vector<Amplitude*> AmpsB = decayInfo->amplitudes_B;

    for(int i=0; i<AmpsA.size(); i++) {


        AmpMap[AmpsA[i]->_uniqueDecayStr] =  std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));

        // printf("Adding Amplitde A:%s\n",AmpsA[i]->_uniqueDecayStr.c_str());

        auto LSvec = AmpsA[i]->_LS;

        for(auto LSIT = LSvec.begin(); LSIT != LSvec.end(); ++LSIT) {
            auto found = std::find_if(components.begin(), components.end(), [LSIT](const PdfBase* L) {
                return (**LSIT)== *(dynamic_cast<const Lineshape*>(L));
            });

            if(found != components.end()) {
                AmpMap[AmpsA[i]->_uniqueDecayStr].first.push_back(std::distance(components.begin(), found));
                // printf("LS %s found at %i\n",(*found)->getName().c_str(),std::distance(components.begin(), found));
            } else {
                components.push_back(*LSIT);
                AmpMap[AmpsA[i]->_uniqueDecayStr].first.push_back(components.size() - 1);
                // printf("Adding LS %s\n",(*LSIT)->getName().c_str());

            }
        }

        auto SFvec = AmpsA[i]->_SF;

        for(auto SFIT = SFvec.begin(); SFIT != SFvec.end(); ++SFIT) {
            auto found = std::find_if(SpinFactors.begin(), SpinFactors.end(), [SFIT](const SpinFactor* S) {
                return (**SFIT) == (*S);
            });

            if(found != SpinFactors.end()) {
                AmpMap[AmpsA[i]->_uniqueDecayStr].second.push_back(std::distance(SpinFactors.begin(), found));
                // printf("SF %s found at %i\n",(*found)->getName().c_str(), std::distance(SpinFactors.begin(), found));
            } else {
                SpinFactors.push_back(*SFIT);
                AmpMap[AmpsA[i]->_uniqueDecayStr].second.push_back(SpinFactors.size() - 1);
                // printf("Adding SF %s\n",(*SFIT)->getName().c_str());

            }
        }


        nPermVec.push_back(AmpsA[i]->_nPerm);
        pindices.push_back(registerParameter(AmpsA[i]->_ar));
        pindices.push_back(registerParameter(AmpsA[i]->_ai));
        ++coeff_counter;
        AmpBuffer.push_back(AmpsA[i]);
        unsigned int flag = 0;
        auto inB = std::find_if(AmpsB.begin(), AmpsB.end(), [AmpsA, i](const Amplitude* A) {
            return *(AmpsA[i]) == (*A);
        });

        if(inB != AmpsB.end()) {
            // printf("Found in AmpsB as well: %s\n", (*inB)->_uniqueDecayStr.c_str());
            flag = 2;
            pindices.push_back(registerParameter((*inB)->_ar));
            pindices.push_back(registerParameter((*inB)->_ai));
            ++coeff_counter;
        }

        ampidxstart.push_back(ampidx.size());
        std::vector<unsigned int> ls = AmpMap[AmpsA[i]->_uniqueDecayStr].first;
        std::vector<unsigned int> sf = AmpMap[AmpsA[i]->_uniqueDecayStr].second;
        ampidx.push_back(ls.size());
        ampidx.push_back(sf.size());
        ampidx.push_back(AmpsA[i]->_nPerm);
        ampidx.push_back(flag);
        ampidx.insert(ampidx.end(), ls.begin(), ls.end());
        ampidx.insert(ampidx.end(), sf.begin(), sf.end());
    }


    for(int i=0; i<AmpsB.size(); i++) {

        unsigned int flag = 1;
        auto inB = std::find_if(AmpBuffer.begin(), AmpBuffer.end(), [AmpsB, i](const Amplitude* A) {
            return *(AmpsB[i]) == (*A);
        });

        if(inB != AmpBuffer.end())
            continue;

        AmpMap[AmpsB[i]->_uniqueDecayStr] =  std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));
        // fprintf("Adding Amplitude B %s\n",AmpsB[i]->_uniqueDecayStr.c_str());

        auto LSvec = AmpsB[i]->_LS;

        for(auto LSIT = LSvec.begin(); LSIT != LSvec.end(); ++LSIT) {
            auto found = std::find_if(components.begin(), components.end(), [LSIT](const PdfBase* L) {
                return (**LSIT)== *(dynamic_cast<const Lineshape*>(L));
            });

            if(found != components.end()) {
                AmpMap[AmpsB[i]->_uniqueDecayStr].first.push_back(std::distance(components.begin(), found));
                // fprintf("LS %s found at %i\n",(*found)->getName().c_str(), std::distance(components.begin(), found));
            } else {
                components.push_back(*LSIT);
                AmpMap[AmpsB[i]->_uniqueDecayStr].first.push_back(components.size() - 1);
                // fprintf("Adding LS %s\n",(*LSIT)->getName().c_str());
            }
        }

        auto SFvec = AmpsB[i]->_SF;

        for(auto SFIT = SFvec.begin(); SFIT != SFvec.end(); ++SFIT) {
            auto found = std::find_if(SpinFactors.begin(), SpinFactors.end(), [SFIT](const SpinFactor* S) {
                return (**SFIT) == (*S);
            });

            if(found != SpinFactors.end()) {
                AmpMap[AmpsB[i]->_uniqueDecayStr].second.push_back(std::distance(SpinFactors.begin(), found));
                // fprintf("SF %s found at %i\n",(*found)->getName().c_str(), std::distance(SpinFactors.begin(), found));
            } else {
                SpinFactors.push_back(*SFIT);
                AmpMap[AmpsB[i]->_uniqueDecayStr].second.push_back(SpinFactors.size() - 1);
                // fprintf("Adding SF %s\n",(*SFIT)->getName().c_str());
            }
        }

        nPermVec.push_back(AmpsB[i]->_nPerm);
        pindices.push_back(registerParameter(AmpsB[i]->_ar));
        pindices.push_back(registerParameter(AmpsB[i]->_ai));
        ++coeff_counter;
        ampidxstart.push_back(ampidx.size());
        std::vector<unsigned int> ls = AmpMap[AmpsB[i]->_uniqueDecayStr].first;
        std::vector<unsigned int> sf = AmpMap[AmpsB[i]->_uniqueDecayStr].second;
        ampidx.push_back(ls.size());
        ampidx.push_back(sf.size());
        ampidx.push_back(AmpsB[i]->_nPerm);
        ampidx.push_back(flag);
        ampidx.insert(ampidx.end(), ls.begin(), ls.end());
        ampidx.insert(ampidx.end(), sf.begin(), sf.end());
    }


    MEMCPY_TO_SYMBOL(AmpIndices, &(ampidxstart[0]), ampidxstart.size()*sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(AmpIndices, &(ampidx[0]), ampidx.size()*sizeof(unsigned int), ampidxstart.size()*sizeof(unsigned int),
                     cudaMemcpyHostToDevice);

    pindices[2] = components.size();
    pindices[3] = SpinFactors.size();
    pindices[4] = AmpMap.size();
    pindices[5] = coeff_counter;

    for(int i = 0; i<components.size(); i++) {
        reinterpret_cast<Lineshape*>(components[i])->setConstantIndex(cIndex);
        pindices.push_back(reinterpret_cast<Lineshape*>(components[i])->getFunctionIndex());
        pindices.push_back(reinterpret_cast<Lineshape*>(components[i])->getParameterIndex());

    }

    for(int i = 0; i < SpinFactors.size(); ++i) {
        pindices.push_back(SpinFactors[i]->getFunctionIndex());
        pindices.push_back(SpinFactors[i]->getParameterIndex());
        SpinFactors[i]->setConstantIndex(cIndex);
    }


    pindices.push_back(efficiency->getFunctionIndex());
    pindices.push_back(efficiency->getParameterIndex());
    components.push_back(efficiency);

    // In case the resolution function needs parameters, this registers them.
    resolution->createParameters(pindices, this);
    GET_FUNCTION_ADDR(ptr_to_TDDP4);
    initialise(pindices);

    Integrator =  new NormIntegrator_TD(parameters);
    redoIntegral = new bool[components.size() - 1];
    cachedMasses = new fptype[components.size() - 1];
    cachedWidths = new fptype[components.size() - 1];

    for(int i = 0; i < components.size() - 1; ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        lscalculators.push_back(new LSCalculator_TD(parameters, i));

    }

    for(int i = 0; i < SpinFactors.size(); ++i) {
        sfcalculators.push_back(new SFCalculator_TD(parameters, i));

    }

    for(int i = 0; i < AmpMap.size(); ++i) {
        AmpCalcs.push_back(new AmpCalc_TD(ampidxstart[i], parameters, nPermVec[i]));
    }

    // fprintf(stderr,"#Amp's %i, #LS %i, #SF %i \n", AmpMap.size(), components.size()-1, SpinFactors.size() );

    std::vector<mcbooster::GReal_t> masses(decayInfo->particle_masses.begin()+1, decayInfo->particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo->particle_masses[0], masses, MCeventsNorm, generation_offset);
    phsp.Generate(mcbooster::Vector4R(decayInfo->particle_masses[0], 0.0, 0.0, 0.0));
    phsp.Unweight();

    auto nAcc = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1 = phsp.GetDaughters(0);
    auto d2 = phsp.GetDaughters(1);
    auto d3 = phsp.GetDaughters(2);
    auto d4 = phsp.GetDaughters(3);

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    auto zip_end = zip_begin + d1.size();
    auto new_end = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    // fprintf("After accept-reject we will keep %.i Events for normalization.\n", (int)nAcc);
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
    VarSet[0] = &norm_M12,
                VarSet[1] = &norm_M34;
    VarSet[2] = &norm_CosTheta12;
    VarSet[3] = &norm_CosTheta34;
    VarSet[4] = &norm_phi;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet);

    norm_SF = mcbooster::RealVector_d(nAcc * SpinFactors.size());
    norm_LS = mcbooster::mc_device_vector<devcomplex<fptype>>(nAcc * (components.size() - 1));
    MCevents = nAcc;


    addSpecialMask(PdfBase::ForceSeparateNorm);
}


// makes the arrays to chache the lineshape values and spinfactors in CachedResSF and the values of the amplitudes in cachedAMPs
// I made the choice to have spinfactors necxt to the values of the lineshape in memory. I waste memory by doing this because a spinfactor is saved as complex
// It would be nice to test if this is better than having the spinfactors stored seperately.
__host__ void TDDP4::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum = 6
    totalEventSize = evtSize;
    assert(totalEventSize >= 3);

    if(cachedResSF)
        delete cachedResSF;

    if(cachedAMPs)
        delete cachedAMPs;

    numEntries = dataSize;
    cachedResSF = new DEVICE_VECTOR<devcomplex<fptype>>(dataSize*(components.size() + SpinFactors.size() -
            1)); //   -1 because 1 component is efficiency
    void* dummy = thrust::raw_pointer_cast(cachedResSF->data());
    MEMCPY_TO_SYMBOL(cResSF_TD, &dummy, sizeof(devcomplex<fptype>*), cacheToUse*sizeof(devcomplex<fptype>*),
                     cudaMemcpyHostToDevice);

    cachedAMPs = new DEVICE_VECTOR<devcomplex<fptype>>(dataSize*(AmpCalcs.size()));
    void* dummy2 = thrust::raw_pointer_cast(cachedAMPs->data());
    MEMCPY_TO_SYMBOL(Amps_TD, &dummy2, sizeof(devcomplex<fptype>*), cacheToUse*sizeof(devcomplex<fptype>*),
                     cudaMemcpyHostToDevice);

    setForceIntegrals();
}

// this is where the actual magic happens. This function does all the calculations!
__host__ fptype TDDP4::normalise() const {
    // fprintf(stderr, "start normalise\n");
    recursiveSetNormalisation(1); // Not going to normalise efficiency,
    // so set normalisation factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);

//check if MINUIT changed any parameters and if so remember that so we know
// we need to recalculate that lineshape and every amp, that uses that lineshape
    for(unsigned int i = 0; i < components.size() - 1; ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(components[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
        components[i]->storeParameters();
    }

    SpinsCalculated = !forceRedoIntegrals;
    forceRedoIntegrals = false;

    //just some thrust iterators for the calculation.
    thrust::constant_iterator<fptype*> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    //Calculate spinfactors only once for normalisation events and real events
    //strided_range is a template implemented in DalitsPlotHelpers.hh
    //it basically goes through the array by increasing the pointer by a certain amount instead of just one step.
    if(!SpinsCalculated) {
        for(int i = 0; i < SpinFactors.size() ; ++i) {
            unsigned int offset = components.size() -1;
            thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                              thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                              strided_range<DEVICE_VECTOR<devcomplex<fptype>>::iterator>(cachedResSF->begin() + offset + i,
                                      cachedResSF->end(),
                                      (components.size() + SpinFactors.size() - 1)).begin(),
                              *(sfcalculators[i]));

            if(!generation_no_norm) {
                thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(), norm_M34.begin(),
                                  norm_CosTheta12.begin(), norm_CosTheta34.begin(), norm_phi.begin()))
                                  , thrust::make_zip_iterator(thrust::make_tuple(norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(),
                                          norm_CosTheta34.end(), norm_phi.end()))
                                  , (norm_SF.begin() + (i * MCevents))
                                  , NormSpinCalculator_TD(parameters, i));
            }
        }

        SpinsCalculated = true;
    }

    // fprintf(stderr, "normalise after spins\n");

    // this calculates the values of the lineshapes and stores them in the array. It is recalculated every time parameters change.
    for(int i = 0; i < components.size() -1 ; ++i) {
        if(redoIntegral[i]) {
            thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                              thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                              strided_range<DEVICE_VECTOR<devcomplex<fptype>>::iterator>(cachedResSF->begin() + i,
                                      cachedResSF->end(),
                                      (components.size() + SpinFactors.size() - 1)).begin(),
                              *(lscalculators[i]));
        }
    }

    // fprintf(stderr, "normalise after LS\n");

    // this is a little messy but it basically checks if the amplitude includes one of the recalculated lineshapes and if so recalculates that amplitude
    std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int>>>::const_iterator AmpMapIt =
        AmpMap.begin();

    for(int i = 0; i < AmpCalcs.size(); ++i) {
        std::vector<unsigned int> redoidx((*AmpMapIt).second.first);
        bool redo = false;

        for(int j = 0; j < redoidx.size(); ++j) {
            if(!redoIntegral[redoidx[j]])
                continue;

            redo = true;
            break;
        }

        if(redo) {
            thrust::transform(eventIndex, eventIndex + numEntries,
                              strided_range<DEVICE_VECTOR<devcomplex<fptype>>::iterator>(cachedAMPs->begin() + i,
                                      cachedAMPs->end(), AmpCalcs.size()).begin(),
                              *(AmpCalcs[i]));
        }
    }

    // fprintf(stderr, "normalise after Amps\n");

    // lineshape value calculation for the normalisation, also recalculated every time parameter change
    if(!generation_no_norm) {
        for(int i = 0; i < components.size() -1 ; ++i) {
            if(!redoIntegral[i])
                continue;

            thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(), norm_M34.begin(),
                              norm_CosTheta12.begin(), norm_CosTheta34.begin(), norm_phi.begin()))
                              , thrust::make_zip_iterator(thrust::make_tuple(norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(),
                                      norm_CosTheta34.end(), norm_phi.end()))
                              , (norm_LS.begin() + (i * MCevents))
                              , NormLSCalculator_TD(parameters, i));
        }
    }

    thrust::constant_iterator<fptype*> normSFaddress(thrust::raw_pointer_cast(norm_SF.data()));
    thrust::constant_iterator<devcomplex<fptype>* > normLSaddress(thrust::raw_pointer_cast(norm_LS.data()));
    thrust::constant_iterator<int> NumNormEvents(MCevents);

    //this does the rest of the integration with the cached lineshape and spinfactor values for the normalization events
    auto ret = 1.0;

    if(!generation_no_norm) {
        thrust::tuple<fptype, fptype, fptype, fptype> dummy(0, 0, 0, 0);
        FourDblTupleAdd MyFourDoubleTupleAdditionFunctor;
        thrust::tuple<fptype, fptype, fptype, fptype> sumIntegral;
        sumIntegral = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, NumNormEvents,
                                               normSFaddress, normLSaddress)),
                                               thrust::make_zip_iterator(thrust::make_tuple(eventIndex + MCevents, NumNormEvents, normSFaddress, normLSaddress)),
                                               *Integrator,
                                               dummy,
                                               MyFourDoubleTupleAdditionFunctor);

        // printf("normalise A2/#evts , B2/#evts: %.5g, %.5g\n",thrust::get<0>(sumIntegral)/MCevents, thrust::get<1>(sumIntegral)/MCevents);
        fptype tau     = host_params[host_indices[parameters + 7]];
        fptype xmixing = host_params[host_indices[parameters + 8]];
        fptype ymixing = host_params[host_indices[parameters + 9]];

        ret = resolution->normalisation(thrust::get<0>(sumIntegral), thrust::get<1>(sumIntegral), thrust::get<2>(sumIntegral),
                                        thrust::get<3>(sumIntegral), tau, xmixing, ymixing);

        //MCevents is the number of normalisation events.
        ret/=MCevents;
    }

    host_normalisation[parameters] = 1.0/ret;
    // printf("end of normalise %f\n", ret);
    return ret;
}

__host__ std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h>
TDDP4::GenerateSig(unsigned int numEvents) {
    std::vector<mcbooster::GReal_t> masses(decayInfo->particle_masses.begin()+1, decayInfo->particle_masses.end());
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
    auto dtime_d             = mcbooster::RealVector_d(numEvents);

    thrust::counting_iterator<unsigned int> index_sequence_begin(0);
    thrust::transform(index_sequence_begin,
                      index_sequence_begin + numEvents,
                      dtime_d.begin(),
                      genUni(genlow, genhigh, generation_offset));


    mcbooster::VariableSet_d VarSet_d(5);
    VarSet_d[0] = &SigGen_M12_d,
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

    auto weights = mcbooster::RealVector_d(phsp.GetWeights());
    phsp.~PhaseSpace();

    auto DS = new mcbooster::RealVector_d(8*numEvents);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 5; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 8);
        thrust::copy(VarSet_d[i]->begin(), VarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 5, DS->end(), 8);
    thrust::copy(eventNumber, eventNumber+numEvents, sr.begin());

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr2(DS->begin() + 6, DS->end(), 8);
    thrust::copy(dtime_d.begin(), dtime_d.end(), sr2.begin());

    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(numEvents, 8);


    generation_no_norm=true; // we need no normalization for generation, but we do need to make sure that norm = 1;
    SigGenSetIndices();
    copyParams();
    normalise();
    setForceIntegrals();
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);

    thrust::device_vector<fptype> results(numEvents);
    thrust::constant_iterator<int> eventSize(8);
    thrust::constant_iterator<fptype*> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);
    thrust::constant_iterator<fptype*> weightAddress(thrust::raw_pointer_cast(weights.data()));
    thrust::constant_iterator<fptype*> dtimeAddress(thrust::raw_pointer_cast(dtime_d.data()));

    MetricTaker evalor(this, getMetricPointer("ptr_to_Prob"));
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEvents, arrayAddress, eventSize)),
                      results.begin(),
                      evalor);
    cudaDeviceSynchronize();

    gooFree(dev_event_array);

    thrust::transform(results.begin(), results.end(), weights.begin(), weights.begin(),
                      thrust::multiplies<mcbooster::GReal_t>());
    mcbooster::BoolVector_d flags(numEvents);

    thrust::counting_iterator<mcbooster::GLong_t> first(0);
    thrust::counting_iterator<mcbooster::GLong_t> last = first + numEvents;

    auto max = thrust::max_element(weights.begin(), weights.end());
    fptype wmax = (fptype)*max;
    thrust::transform(first, last, weights.begin(), flags.begin(), mcbooster::FlagAcceptReject(wmax, generation_offset));

    //   printf("Offset: %i und wmax:%.5g\n",generation_offset, wmax );
    // for (int i = 0; i < dtime_d.size(); ++i)
    // {
    // printf("%i, %s, %.5g, %.5g, %.5g, %.5g, %.5g, %.5g, %.5g\n",i, (bool)flags[i] ? "true" : "false", (double)weights[i],  (double)dtime_d[i], (double)SigGen_M12_d[i], (double)SigGen_M34_d[i], (double)SigGen_CosTheta12_d[i], (double)SigGen_CosTheta34_d[i], (double)SigGen_phi_d[i]);
    // }

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h = mcbooster::BoolVector_h(flags);
    cudaDeviceSynchronize();

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}

SFCalculator_TD::SFCalculator_TD(int pIdx, unsigned int sf_idx)
    : _spinfactor_i(sf_idx)
    , _parameters(pIdx)
{}

__device__ devcomplex<fptype> SFCalculator_TD::operator()(thrust::tuple<int, fptype*, int> t) const {

    int evtNum = thrust::get<0>(t);
    fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
    int parameter_i = 12 + (2 * indices[6]) + (indices[3] * 2) + (_spinfactor_i * 2)
                      ; // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i+1];

    fptype m12 = evt[indices[2 + indices[0]]];
    fptype m34 = evt[indices[3 + indices[0]]];
    fptype cos12 = evt[indices[4 + indices[0]]];
    fptype cos34 = evt[indices[5 + indices[0]]];
    fptype phi = evt[indices[6 + indices[0]]];

    fptype vecs[16];
    get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);
    // printf("%i, %i, %f, %f, %f, %f, %f \n",evtNum, thrust::get<2>(t), m12, m34, cos12, cos34, phi );
    // printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
    // printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
    // printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
    // printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);

    spin_function_ptr func = reinterpret_cast<spin_function_ptr>(device_function_table[functn_i]);
    fptype sf = (*func)(vecs, paramIndices+params_i);
    // printf("SpinFactors %i : %.7g\n",_spinfactor_i, sf );
    return devcomplex<fptype>(sf, 0);
}

NormSpinCalculator_TD::NormSpinCalculator_TD(int pIdx, unsigned int sf_idx)
    : _spinfactor_i(sf_idx)
    , _parameters(pIdx)
{}

__device__ fptype NormSpinCalculator_TD::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
const {

    unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
    int parameter_i = 12 + (2 * indices[6]) + (indices[3] * 2) + (_spinfactor_i * 2)
                      ; // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i+1];

    fptype m12    = (thrust::get<0>(t));
    fptype m34    = (thrust::get<1>(t));
    fptype cos12  = (thrust::get<2>(t));
    fptype cos34  = (thrust::get<3>(t));
    fptype phi    = (thrust::get<4>(t));

    fptype vecs[16];
    get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);

//   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,0, vecs[0], vecs[1], vecs[2], vecs[3]);
//   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,1, vecs[4], vecs[5], vecs[6], vecs[7]);
//   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,2, vecs[8], vecs[9], vecs[10], vecs[11]);
//   printf("evt %i vec%i %.5g, %.5g, %.5g, %.5g\n", evtNum,3, vecs[12], vecs[13], vecs[14], vecs[15]);
// // }
    spin_function_ptr func = reinterpret_cast<spin_function_ptr>(device_function_table[functn_i]);
    fptype sf = (*func)(vecs, paramIndices+params_i);

    // printf("NormSF evt:%.5g, %.5g, %.5g, %.5g, %.5g\n", m12, m34, cos12, cos34, phi);
    // printf("NormSF %i, %.7g\n",_spinfactor_i, sf );
    // THREAD_SYNCH
    return sf;
}



LSCalculator_TD::LSCalculator_TD(int pIdx, unsigned int res_idx)
    : _resonance_i(res_idx)
    , _parameters(pIdx)
{}

__device__ devcomplex<fptype> LSCalculator_TD::operator()(thrust::tuple<int, fptype*, int> t) const {
    // Calculates the BW values for a specific resonance.
    devcomplex<fptype> ret;

    int evtNum = thrust::get<0>(t);
    fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
    int parameter_i = 12 + (2 * indices[6]) + (_resonance_i * 2)
                      ; // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i+1];
    unsigned int pair = (paramIndices+params_i)[5];

    fptype m1  = functorConstants[indices[1] + 2];
    fptype m2  = functorConstants[indices[1] + 3];
    fptype m3  = functorConstants[indices[1] + 4];
    fptype m4  = functorConstants[indices[1] + 5];

    fptype m12 = evt[indices[2 + indices[0]]];
    fptype m34 = evt[indices[3 + indices[0]]];
    fptype cos12 = evt[indices[4 + indices[0]]];
    fptype cos34 = evt[indices[5 + indices[0]]];
    fptype phi = evt[indices[6 + indices[0]]];

    if(pair < 2) {
        fptype mres = pair==0 ? m12 : m34;
        fptype d1 = pair==0 ? m1 : m3;
        fptype d2 = pair==0 ? m2 : m4;
        ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
        // printf("LS_nt %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    } else {
        fptype vecs[16];
        get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
        // printf("LS %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );

    }

    // printf("LS %i: %.7g, %.7g, %.7g, %.7g, %.7g \n",evtNum, m12, m34, cos12, cos34, phi );

    //if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
    //printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    //printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    // printf("BW_%i : %f %f\n", _resonance_i, ret.real, ret.imag);

    return ret;
}

NormLSCalculator_TD::NormLSCalculator_TD(int pIdx, unsigned int res_idx)
    : _resonance_i(res_idx)
    , _parameters(pIdx)
{}

__device__ devcomplex<fptype> NormLSCalculator_TD::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
const {
    // Calculates the BW values for a specific resonance.
    devcomplex<fptype> ret;

    unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
    int parameter_i = 12 + (2 * indices[6]) + (_resonance_i *
                      2); // Find position of this resonance relative to DALITZPLOT start
    unsigned int functn_i = indices[parameter_i];
    unsigned int params_i = indices[parameter_i+1];
    unsigned int pair = (paramIndices+params_i)[5];

    fptype m1  = functorConstants[indices[1] + 2];
    fptype m2  = functorConstants[indices[1] + 3];
    fptype m3  = functorConstants[indices[1] + 4];
    fptype m4  = functorConstants[indices[1] + 5];

    fptype m12    = (thrust::get<0>(t));
    fptype m34    = (thrust::get<1>(t));
    fptype cos12  = (thrust::get<2>(t));
    fptype cos34  = (thrust::get<3>(t));
    fptype phi    = (thrust::get<4>(t));


    if(pair < 2) {
        fptype mres = pair==0 ? m12 : m34;
        fptype d1 = pair==0 ? m1 : m3;
        fptype d2 = pair==0 ? m2 : m4;
        ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
    } else {
        fptype vecs[16];
        get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
    }

    // printf("NormLS %f, %f, %f, %f, %f \n",m12, m34, cos12, cos34, phi );
    // printf("%i, %i, %i, %i, %i \n",numLS, numSF, numAmps, offset, evtNum );
    // printf("NLS %i, %f, %f\n",_resonance_i,ret.real, ret.imag);

    //printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    //printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    // THREAD_SYNCH
    return ret;
}

AmpCalc_TD::AmpCalc_TD(unsigned int AmpIdx, unsigned int pIdx, unsigned int nPerm)
    : _nPerm(nPerm)
    , _AmpIdx(AmpIdx)
    , _parameters(pIdx)
{}


__device__ devcomplex<fptype> AmpCalc_TD::operator()(thrust::tuple<int, fptype*, int> t) const {
    unsigned int* indices = paramIndices + _parameters;
    unsigned int cacheToUse = indices[2];
    unsigned int totalLS = indices[3];
    unsigned int totalSF = indices[4];
    unsigned int totalAMP = indices[5];
    unsigned int offset = totalLS + totalSF;
    unsigned int numLS = AmpIndices[totalAMP + _AmpIdx];
    unsigned int numSF = AmpIndices[totalAMP + _AmpIdx + 1];
    unsigned int evtNum = thrust::get<0>(t);



    devcomplex<fptype> returnVal(0, 0);
    unsigned int SF_step = numSF/_nPerm;
    unsigned int LS_step = numLS/_nPerm;

    for(int i = 0; i < _nPerm; ++i) {
        devcomplex<fptype> ret(1, 0);

        for(int j = i*LS_step; j < (i+1)*LS_step; ++j) {
            ret *= (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 + j]]);
            // printf("Lineshape %i = (%.7g, %.7g)\n", j, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 + j]]).real, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 + j]]).imag);

        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i*SF_step; j < (i+1)*SF_step; ++j) {
            ret *= (cResSF_TD[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 4 + numLS + j]].real);
            // printf(" SF = %.7g\n", (cResSF_TD[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 4 + numLS + j]].real));

        }

        // printf("Lineshape Product * SF = (%.7g, %.7g)\n", ret.real, ret.imag);

        returnVal += ret;
    }

    returnVal *= (1/SQRT((fptype)(_nPerm)));
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real, returnVal.imag);
    return  returnVal;
}


NormIntegrator_TD::NormIntegrator_TD(unsigned int pIdx)
    : _parameters(pIdx)
{}


__device__ thrust::tuple<fptype, fptype, fptype, fptype> NormIntegrator_TD::operator()(
    thrust::tuple<int, int, fptype*, devcomplex<fptype>*> t) const {
    unsigned int* indices = paramIndices + _parameters;
    unsigned int totalAMP = indices[5];

    unsigned int evtNum = thrust::get<0>(t);
    unsigned int MCevents = thrust::get<1>(t);
    fptype* SFnorm = thrust::get<2>(t) + evtNum;
    devcomplex<fptype>* LSnorm = thrust::get<3>(t) + evtNum;

    devcomplex<fptype> AmpA(0, 0);
    devcomplex<fptype> AmpB(0, 0);
    fptype amp_real_A, amp_imag_A, amp_real_B, amp_imag_B;

    int k = 0;

    for(int amp = 0; amp < totalAMP; ++amp) {
        unsigned int ampidx =  AmpIndices[amp];
        unsigned int numLS = AmpIndices[totalAMP + ampidx];
        unsigned int numSF = AmpIndices[totalAMP + ampidx + 1];
        unsigned int nPerm = AmpIndices[totalAMP + ampidx + 2];
        unsigned int flag = AmpIndices[totalAMP + ampidx + 3];
        unsigned int SF_step = numSF/nPerm;
        unsigned int LS_step = numLS/nPerm;
        devcomplex<fptype> ret2(0, 0);
        // printf("%i, %i, %i, %i, %i, %i, %i, %i, %i, %f\n",ampidx, amp, numLS, numSF, nPerm,AmpIndices[totalAMP + ampidx + 4 + 0], AmpIndices[totalAMP + ampidx + 4 + 1], AmpIndices[totalAMP + ampidx + 4 + 2], AmpIndices[totalAMP + ampidx + 4 + 3], (1/SQRT((fptype)(nPerm))) );

        for(int j = 0; j < nPerm; ++j) {
            devcomplex<fptype> ret(1, 0);

            for(int i = j*LS_step; i < (j+1)*LS_step; ++i) {
                devcomplex<fptype> matrixelement(LSnorm[AmpIndices[totalAMP + ampidx + 4 + i] * MCevents]);
                // printf("Norm BW %i, %.5g, %.5g\n",AmpIndices[totalAMP + ampidx + 4 + i] , matrixelement.real, matrixelement.imag);
                ret *= matrixelement;

            }

            for(int i = j*SF_step; i < (j+1)*SF_step; ++i) {
                fptype matrixelement = (SFnorm[AmpIndices[totalAMP + ampidx + 4 + numLS + i] * MCevents]);
                // printf("Norm SF %i, %.5g\n",AmpIndices[totalAMP + ampidx + 4 + i] , matrixelement);
                ret *= matrixelement;

            }

            ret2 += ret;
        }

        ret2 *= (1/SQRT((fptype)(nPerm)));
        // printf("Result Amplitude %i, %i, %.5g, %.5g\n",flag, amp, ret2.real, ret2.imag);

        switch(flag) {
        case 0:
            amp_real_A = cudaArray[indices[12 + 2*(amp+k)]];
            amp_imag_A = cudaArray[indices[13 + 2*(amp+k)]];
            // printf("DEV0: %.5g, %.5g, %.5g +i %.5g \n", amp_real_A, amp_imag_A, ret2.real, ret2.imag);
            AmpA += ((ret2).multiply(amp_real_A, amp_imag_A));
            break;

        case 1:
            amp_real_B = cudaArray[indices[12 + 2*(amp+k)]];
            amp_imag_B = cudaArray[indices[13 + 2*(amp+k)]];
            // printf("DEV1: %.5g, %.5g, %.5g +i %.5g \n", amp_real_B, amp_imag_B, ret2.real, ret2.imag);
            AmpB += ((ret2).multiply(amp_real_B, amp_imag_B));
            break;

        case 2:
            amp_real_A = cudaArray[indices[12 + 2*(amp+k)]];
            amp_imag_A = cudaArray[indices[13 + 2*(amp+k)]];
            // printf("DEV2.1: %.5g, %.5g, %.5g +i %.5g \n", amp_real_A, amp_imag_A, ret2.real, ret2.imag);
            AmpA += ((ret2).multiply(amp_real_A, amp_imag_A));
            ++k;
            amp_real_B = cudaArray[indices[12 + 2*(amp+k)]];
            amp_imag_B = cudaArray[indices[13 + 2*(amp+k)]];
            // printf("DEV2.2: %.5g, %.5g, %.5g +i %.5g \n", amp_real_B, amp_imag_B, ret2.real, ret2.imag);
            AmpB += ((ret2).multiply(amp_real_B, amp_imag_B));
            break;
        }
    }

    fptype _SqWStoRSrate = cudaArray[indices[10]];
    AmpA *= _SqWStoRSrate;

    auto AmpAB = AmpA*conj(AmpB);
    // printf("%.5g %.5g %.5g %.5g\n", norm2(AmpA), norm2(AmpB), AmpAB.real, AmpAB.imag);
    return thrust::tuple<fptype, fptype, fptype, fptype>(norm2(AmpA), norm2(AmpB), AmpAB.real, AmpAB.imag);
}
