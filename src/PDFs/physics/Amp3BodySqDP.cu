#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3BodySqDP.h>
#include <goofit/PDFs/physics/Amp3BodyBase.h>
#include <goofit/PDFs/physics/detail/Dim2.h>
#include <goofit/PDFs/physics/detail/SpecialSqDpResonanceCalculator.h>
#include <goofit/PDFs/physics/detail/SpecialSqDpResonanceIntegrator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/detail/Complex.h>

#include <thrust/copy.h>
#include <thrust/transform_reduce.h>
#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/system/omp/execution_policy.h>
#include <thrust/system/cpp/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>
#include <array>
#include <vector>



namespace GooFit {


__host__ __device__   auto inSqDalitz(const fptype &mprime,const fptype &thetaprime) -> bool{
    return (mprime>0.0 && mprime<1.0)&&(thetaprime>0.0 && thetaprime<1.0);
}

__host__ __device__  auto calc_mprime(const fptype &m12, const fptype &m_mother, const fptype &m1, const fptype &m2, const fptype &m3)->fptype{
    fptype min = m1+m2;
    fptype max = m_mother-m3;
    fptype mprime = (2*(m12 - min)/(max-min)) - 1.0;

    if(mprime<-1.)
        mprime=-1.;
    
    if(mprime>1.)
        mprime=1.;

    return acos(mprime)/M_PI;
}

__host__ __device__  auto calc_thetaprime(const fptype &m12,const fptype &m13, const fptype &m_mother, const fptype &m1, const fptype &m2, const fptype &m3)->fptype{
    fptype m13Sq = m13*m13;
    fptype m12Sq = m12*m12;
    fptype m_motherSq = m_mother*m_mother;
    fptype m1Sq = m1*m1;
    fptype m2Sq = m2*m2;
    fptype m3Sq = m3*m3;

    fptype EiCmsij = (m12Sq - m2Sq + m1Sq)/(2.0*m12);
    fptype EkCmsij = (m_motherSq - m12Sq - m3Sq)/(2.0*m12);

    //printf("EiCmsij=%f \t EkCmsij=%f \n",m12, m13);

    fptype qi = EiCmsij*EiCmsij - m1Sq;
    qi = qi>0. ? sqrt(qi) : 0.;

    fptype qk = EkCmsij*EkCmsij - m3Sq;
    qk = qk>0. ?  sqrt(qk) : 0.;
    
    fptype coshel = (m13Sq - m1Sq - m3Sq - 2.0*EiCmsij*EkCmsij)/(2.0*qi*qk);

    if(coshel<-1.)
        coshel=-1.;
    
    if(coshel>1.)
        coshel=1.;

    fptype thetaprime = acos(coshel)/M_PI;

    // if(thetaprime>0.5)
    //     thetaprime=1.0-thetaprime;

    return thetaprime;
}

__host__ __device__   auto calc_m12(const fptype &mprime, const fptype &m_mother, const fptype &m1, const fptype &m2, const fptype &m3)->fptype{

    fptype m12 = 0.5*( (m_mother-m3) - (m1+m2) )*(1.0 + cos(M_PI*mprime)) + (m1+m2);

    return m12;
}

__host__ __device__   auto calc_m13(const fptype &m12, const fptype &cos_12, const fptype &m_mother, const fptype &m1, const fptype &m2, const fptype &m3)->fptype{
    
    fptype m12Sq = m12*m12;
    fptype m_motherSq = m_mother*m_mother;
    fptype m1Sq = m1*m1;
    fptype m2Sq = m2*m2;
    fptype m3Sq = m3*m3;

    fptype EiCmsij = (m12Sq - m2Sq + m1Sq)/(2.0*m12);
    fptype EkCmsij = (m_motherSq - m12Sq - m3Sq)/(2.0*m12);

    fptype qi = EiCmsij*EiCmsij - m1Sq;
    qi = qi>0. ? sqrt(qi) : 0.;

    fptype qk = EkCmsij*EkCmsij - m3Sq;
    qk = qk>0. ? sqrt(qk)  : 0.;

    //printf("coshel = %.2f \n",cos_12);
    
    fptype m13Sq = m1Sq + m3Sq + 2.0*EiCmsij*EkCmsij - 2.0*qi*qk*cos_12;

    return sqrt(m13Sq);
}


__host__ __device__  auto calc_SqDp_Jacobian(const fptype &mprime ,const fptype &thetaprime, const fptype &m_mother, const fptype &m1, const fptype &m2, const fptype &m3)->fptype{

    fptype m12 = calc_m12(mprime,m_mother,m1,m2,m3);
    fptype m12Sq = m12*m12;

    fptype m_motherSq = m_mother*m_mother;
    fptype m1Sq = m1*m1;
    fptype m2Sq = m2*m2;
    fptype m3Sq = m3*m3;

    fptype EiCmsij = (m12Sq - m2Sq + m1Sq)/(2.0*m12);
    fptype EkCmsij = (m_motherSq - m12Sq - m3Sq)/(2.0*m12);
    
    fptype qi = EiCmsij*EiCmsij - m1Sq;
    qi = qi>0. ? sqrt(qi) : 0.;

    fptype qk = EkCmsij*EkCmsij - m3Sq;
     qk = qk>0. ? sqrt(qk)  : 0.;
    
    fptype deriv1 = 0.5*M_PI*((m_mother-m3) - (m1+m2))*sin(M_PI*mprime);
    fptype deriv2 = M_PI*sin(M_PI*thetaprime);

    fptype jacobian = 4.0*qi*qk*m12*deriv1*deriv2;

    return jacobian;
}

struct prg
{
    fptype a, b;

    __host__ __device__
    prg(fptype _a=0., fptype _b=1.) : a(_a), b(_b) {};

    __host__ __device__
        fptype operator()(const unsigned int n) const
        {
            thrust::default_random_engine rng(n);
            thrust::uniform_real_distribution<fptype> dist(a, b);
            rng.discard(n);

            return dist(rng);
        }
};



__host__  void genNormFakeEvents(size_t n){
   
    thrust::device_vector<fptype> mprime(n);
    thrust::device_vector<fptype> thprime(n);

    thrust::counting_iterator<unsigned int> index_sequence(0);    

    thrust::transform(index_sequence,
            index_sequence + n,
            mprime.begin(),
            prg(0.,1.));

    thrust::counting_iterator<unsigned int> index_sequence2(324234234);
    thrust::transform(index_sequence2,
            index_sequence2 + n,
            thprime.begin(),
            prg(0.,1.));

    mcbooster::VariableSet_d GooVarSet_d(2);
    GooVarSet_d[0] = &mprime;
    GooVarSet_d[1] = &thprime;

    auto DS = new mcbooster::RealVector_d(3 * n);

    thrust::counting_iterator<unsigned int> eventNumber(0);

#pragma unroll
    for(int i = 0; i < 2; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 3);
        thrust::copy(GooVarSet_d[i]->begin(), GooVarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 2, DS->end(), 3);
    thrust::copy(eventNumber, eventNumber + n, sr.begin());

    dev_fake_event_array = thrust::raw_pointer_cast(DS->data());
   

    
}



// Functor used for fit fraction sum
struct CoefSumFunctor {
    fpcomplex coef_i;
    fpcomplex coef_j;

    CoefSumFunctor(fpcomplex coef_i, fpcomplex coef_j)
        : coef_i(coef_i)
        , coef_j(coef_j) {}

    __device__ auto operator()(thrust::tuple<fpcomplex, fpcomplex> val) -> fptype {
        return (coef_i * thrust::conj<fptype>(coef_j) * thrust::get<0>(val) * thrust::conj<fptype>(thrust::get<1>(val)))
            .real();
    }
};

constexpr int resonanceOffset_DP = 4; // Offset of the first resonance into the parameter index array
// Offset is number of parameters, constant index, number of resonances (not calculable
// from nP because we don't know what the efficiency might need), and cache index. Efficiency
// parameters are after the resonance information.

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!

// NOTE: This is does not support ten instances (ten threads) of resoncances now, only one set of resonances.
// this needs to be large enough to hold all samples
__device__ fpcomplex *cSqDpResonances[16 * 20];
fptype cResonancesIntegrals[20][16];

__device__ inline auto parIndexFromResIndex_DP(int resIndex) -> int {
    return resonanceOffset_DP + resIndex * resonanceSize;
}


__device__ auto device_SqDalitzPlot(fptype *evt, ParameterContainer &pc) -> fptype {
    int num_obs = pc.getNumObservables();
    int id_mprime  = pc.getObservable(0);
    int id_thetaprime  = pc.getObservable(1);
    int id_num  = pc.getObservable(2);

    fptype mprime = RO_CACHE(evt[id_mprime]);
    fptype thetaprime = RO_CACHE(evt[id_thetaprime]);

    unsigned int numResonances = pc.getConstant(0);
    unsigned int cacheToUse    = pc.getConstant(1);

    if(!inSqDalitz(mprime, thetaprime)) {
        pc.incrementIndex(1, numResonances * 2, 2, num_obs, 1);

        // loop over resonances and efficiency functions
        for(int i = 0; i < numResonances; i++)
            pc.incrementIndex();

        // increment the efficiency function
        pc.incrementIndex();
        return 0;
    }

    fptype evtIndex = RO_CACHE(evt[id_num]);

    auto evtNum = static_cast<int>(floor(0.5 + evtIndex));

    fpcomplex totalAmp(0, 0);

    for(int i = 0; i < numResonances; ++i) {
        for(int j = 0; j < numResonances; ++j) {
            fpcomplex amp_i = fpcomplex(pc.getParameter(i * 2), pc.getParameter(i * 2 + 1));
            fpcomplex me_i = RO_CACHE(cSqDpResonances[i + (16 * cacheToUse)][evtNum]);
            fpcomplex amp_j = fpcomplex(pc.getParameter(j * 2), -pc.getParameter(j * 2 + 1));
            fpcomplex me_j = RO_CACHE(cSqDpResonances[j + (16 * cacheToUse)][evtNum]);
            totalAmp += amp_i*amp_j*me_i*conj(me_j);
        }
    }

    fptype ret = thrust::abs(totalAmp);


    // for(int i = 0; i < numResonances; ++i) {
    //         fpcomplex amp_i = fpcomplex(pc.getParameter(i * 2), pc.getParameter(i * 2 + 1));
    //         fpcomplex me_i = RO_CACHE(cSqDpResonances[i + (16 * cacheToUse)][evtNum]);
    //         totalAmp += amp_i*me_i;
    // }

    // fptype ret = thrust::norm(totalAmp);
   
    pc.incrementIndex(1, numResonances * 2, 2, num_obs, 1);

    // loop to efficiency idx
    for(int i = 0; i < numResonances; i++)
        pc.incrementIndex();

    fptype eff = callFunction(evt, pc);
  
    //fptype jacobian = calc_SqDp_Jacobian(mprime, thetaprime, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass);
    ret *= eff;

    // printf("likelihood=%f eff=%f \n",ret, eff);

    return ret;
}

int Amp3BodySqDP::cacheCount                         = 0;
__device__ device_function_ptr ptr_to_SqDalitzPlot = device_SqDalitzPlot;

__host__ Amp3BodySqDP::Amp3BodySqDP(
    std::string n, Observable mprime, Observable thetaprime, EventNumber eventNumber, DecayInfo3 decay, GooPdf *efficiency)
    : Amp3BodyBase("Amp3BodySqDP", n, mprime, thetaprime, eventNumber)
    , decayInfo(decay)
    , _mprime(mprime)
    , _thetaprime(thetaprime)
    , _eventNumber(eventNumber)
    , dalitzNormRange(nullptr)
    //, cachedWaves(0)
    , integrals(nullptr)
    , integrals_ff(nullptr)
    , NumNormEvents(decay.SetSizeNormSample)
    , forceRedoIntegrals(true)
    , totalEventSize(3) // Default 3 = mprime, thetaprime, evtNum
    , cacheToUse(0)
    , integrators(nullptr)
    , integrators_ff(nullptr)
    , calculators(nullptr) {
    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    // Passing values to the defined constants.  Rather than push into list, which means each resonance
    MEMCPY_TO_SYMBOL(c_motherMass, &decay.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decay.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decay.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decay.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decay.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_mother_meson_radius, &decay.mother_meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    

    // registered to 0 position
    registerConstant(decayInfo.resonances.size());

    // TODO increase after registerConstant?
    cacheToUse = cacheCount++;
    // registered to 1 position
    registerConstant(cacheToUse);

    for(auto &resonance : decayInfo.resonances) {
        // registering 2 parameters
        registerParameter(resonance->amp_real);
        registerParameter(resonance->amp_imag);
        components.push_back(resonance);
    }

    components.push_back(efficiency);

    registerFunction("ptr_to_SqDalitzPlot", ptr_to_SqDalitzPlot);

    initialize();

    redoIntegral = new bool[decayInfo.resonances.size()];
    cachedMasses = new fptype[decayInfo.resonances.size()];
    cachedWidths = new fptype[decayInfo.resonances.size()];
    integrals    = new fpcomplex **[decayInfo.resonances.size()];
    integrators  = new SpecialSqDpResonanceIntegrator **[decayInfo.resonances.size()];
    integrals_ff    = new fpcomplex **[decayInfo.resonances.size()];
    integrators_ff  = new SpecialSqDpResonanceIntegrator **[decayInfo.resonances.size()];
    calculators  = new SpecialSqDpResonanceCalculator *[decayInfo.resonances.size()];

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        integrators[i]  = new SpecialSqDpResonanceIntegrator *[decayInfo.resonances.size()];
        calculators[i]  = new SpecialSqDpResonanceCalculator(parameters, i);
        integrals[i]    = new fpcomplex *[decayInfo.resonances.size()];

        integrals_ff[i]    = new fpcomplex *[decayInfo.resonances.size()];
        integrators_ff[i]  = new SpecialSqDpResonanceIntegrator *[decayInfo.resonances.size()];

        for(int j = 0; j < decayInfo.resonances.size(); ++j) {
            integrals[i][j]   = new fpcomplex(0, 0);
            integrators[i][j] = new SpecialSqDpResonanceIntegrator(parameters, i, j);
            integrals_ff[i][j]   = new fpcomplex(0, 0);
            integrators_ff[i][j] = new SpecialSqDpResonanceIntegrator(parameters, i, j);
        }
    }

 

    
    setSeparateNorm();
    
}

void Amp3BodySqDP::populateArrays() {
    PdfBase::populateArrays();

    // save our efficiency function.  Resonance's are saved first, then the efficiency function.  Take -1 as efficiency!
    efficiencyFunction = host_function_table.size() - 1;
}
__host__ void Amp3BodySqDP::setDataSize(unsigned int dataSize, unsigned int evtSize, unsigned int offset) {
    // Default 3 is mprime, thetaprime, evtNum
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} must be 3 or more", totalEventSize);

    // if (cachedWaves) delete cachedWaves;
    if(cachedWaves[0]) {
        for(auto &cachedWave : cachedWaves) {
            delete cachedWave;
            cachedWave = nullptr;
        }
    }

    numEntries  = dataSize;
    eventOffset = offset;

    for(int i = 0; i < 16; i++) {
#ifdef GOOFIT_MPI
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(m_iEventsPerTask);
#else
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(dataSize);
#endif
        void *dummy = thrust::raw_pointer_cast(cachedWaves[i]->data());
        MEMCPY_TO_SYMBOL(cSqDpResonances,
                         &dummy,
                         sizeof(fpcomplex *),
                         ((16 * cacheToUse) + i) * sizeof(fpcomplex *),
                         cudaMemcpyHostToDevice);
    }
  

    setForceIntegrals();

    if(dev_fake_event_array==nullptr){
        std::cout << "START: Fill normalization data array with " << NumNormEvents << " events. \n";
        genNormFakeEvents(NumNormEvents);
        cudaDeviceSynchronize();
        std::cout << "END: Done! \n";
    }

}

__host__ auto Amp3BodySqDP::normalize() -> fptype {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    // we need to update the normal here, as values are used at this point.
    host_normalizations.sync(d_normalizations);
    
    int totalBins = _mprime.getNumBins() * _thetaprime.getNumBins();
    double binSizeFactor = 1;
    binSizeFactor *= _mprime.getBinSize();
    binSizeFactor *= _thetaprime.getBinSize();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));
    }

    // This line runs once
    static std::array<fptype, 6> host_norms{{0, 0, 0, 0, 0, 0}};

    std::array<fptype, 6> current_host_norms{{_mprime.getLowerLimit(),
                                              _mprime.getUpperLimit(),
                                              static_cast<fptype>(_mprime.getNumBins()),
                                              _thetaprime.getLowerLimit(),
                                              _thetaprime.getUpperLimit(),
                                              static_cast<fptype>(_thetaprime.getNumBins())}};

    if(host_norms != current_host_norms) {
        host_norms = current_host_norms;
    }
    MEMCPY(dalitzNormRange, host_norms.data(), 6 * sizeof(fptype), cudaMemcpyHostToDevice);
    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(decayInfo.resonances[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    forceRedoIntegrals = false;

    thrust::constant_iterator<fptype *> arrayAddress(dalitzNormRange);
    thrust::counting_iterator<int> binIndex(0);
    thrust::counting_iterator<int> fakeEvtIndex(0);
   
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<fptype *> fakedataArray(dev_fake_event_array);

    thrust::constant_iterator<int> fakeeventSize(totalEventSize);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(eventOffset);

    //Calculate Resonances integrals without multiply by eff and calculate cache
    //These intergrals are used to normalize the RBWs
    size_t n_res=decayInfo.resonances.size();
    thrust::device_vector<fptype> device_integrals(n_res);
    for(int i = 0; i < n_res; ++i) {
                if((!redoIntegral[i]))
                    continue;

                integrators[i][i]->setDalitzIndex(getFunctionIndex());
                integrators[i][i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
                integrators[i][i]->setEfficiencyIndex(decayInfo.resonances[i]->getFunctionIndex());
                integrators[i][i]->setNoEff();
                thrust::constant_iterator<int> effFunc(efficiencyFunction);
                fpcomplex dummy(0, 0);
                thrust::plus<fpcomplex> complexSum;
                (*(integrals[i][i])) = thrust::transform_reduce(
                    thrust::make_zip_iterator(thrust::make_tuple(fakeEvtIndex, fakedataArray, fakeeventSize, effFunc)),
                    thrust::make_zip_iterator(thrust::make_tuple(fakeEvtIndex + NumNormEvents, arrayAddress, fakeeventSize, effFunc)),
                    *(integrators[i][i]),
                    dummy,
                    complexSum); 

                //Check if integral is 0.
                fptype integral = ((*(integrals[i][i])).real()/NumNormEvents);
                
                if(integral<=0.)
                    GOOFIT_ERROR("ERROR: Integral %d = %f.2f", i, integral); 
              
                fptype norm = 1./integral;

                cResonancesIntegrals[cacheToUse][i] = norm;
                //printf("Res  %d -  norm = %f \n",i,1./norm);

                calculators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
                calculators[i]->setDalitzIndex(getFunctionIndex());
            
                thrust::transform(
                    thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                    thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                    strided_range<thrust::device_vector<fpcomplex>::iterator>(cachedWaves[i]->begin(), cachedWaves[i]->end(), 1).begin(),
                    *(calculators[i])
                );

                
                thrust::transform(
                    cachedWaves[i]->begin(),
                    cachedWaves[i]->end(), 
                    thrust::make_constant_iterator(sqrt(norm)),
                    cachedWaves[i]->begin(),
                    thrust::multiplies<fpcomplex>());
               
                  
    }
    
    //Calculate All Integrals with Eff
    for(int i = 0; i < n_res; ++i) {
        for(int j = 0; j < n_res; ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j]))
                continue;
        
            integrators[i][j]->setDalitzIndex(getFunctionIndex());
            integrators[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            integrators[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());
            thrust::constant_iterator<int> effFunc(efficiencyFunction);
            fpcomplex dummy(0, 0);
            thrust::plus<fpcomplex> complexSum;
            (*(integrals[i][j])) = thrust::transform_reduce(
                thrust::system::cuda::par,
                thrust::make_zip_iterator(thrust::make_tuple(fakeEvtIndex, fakedataArray, fakeeventSize, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(fakeEvtIndex + NumNormEvents, arrayAddress, fakeeventSize, effFunc)),
                *(integrators[i][j]),
                dummy,
                complexSum);
        }
    }

    fpcomplex sumIntegral(0, 0);
    for(unsigned int i = 0; i < n_res; ++i) {
        fpcomplex amplitude_i(host_parameters[parametersIdx + i * 2 + 1], host_parameters[parametersIdx + i * 2 + 2]);
        for(unsigned int j = 0; j < n_res ; ++j) {
            fpcomplex amplitude_j(host_parameters[parametersIdx + j * 2 + 1],
                                  -host_parameters[parametersIdx + j * 2 + 2]);

            fptype fNorm_i = cResonancesIntegrals[cacheToUse][i];
            fptype fNorm_j = cResonancesIntegrals[cacheToUse][j];

            sumIntegral += amplitude_i * amplitude_j * ((*(integrals[i][j]))/NumNormEvents)*sqrt(fNorm_i)*sqrt(fNorm_j);
            
       
        }
    }

    fptype ret           = sumIntegral.real(); // That complex number is a square, so it's fully real
  
    //ret /= totalFakeEvent;
    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    return ret;
}

__host__ auto Amp3BodySqDP::sumCachedWave(size_t i) const -> fpcomplex {
    const thrust::device_vector<fpcomplex> &vec = getCachedWaveNoCopy(i);

    fpcomplex ret = thrust::reduce(vec.begin(), vec.end(), fpcomplex(0, 0), thrust::plus<fpcomplex>());

    return ret;
}

__host__ auto Amp3BodySqDP::getCachedWave(size_t i) const -> const std::vector<std::complex<fptype>> {
    // TODO: This calls itself immediately ?
    auto ret_thrust = getCachedWave(i);
    std::vector<std::complex<fptype>> ret(ret_thrust.size());
    thrust::copy(ret_thrust.begin(), ret_thrust.end(), ret.begin());
    return ret;
}

__host__ auto Amp3BodySqDP::fit_fractions(bool print) -> std::vector<std::vector<fptype>> {
     recursiveSetNormalization(1.0);

    host_normalizations.sync(d_normalizations);
    
    size_t n_res     = getDecayInfo().resonances.size();
    size_t totalBins = _mprime.getNumBins() * _thetaprime.getNumBins();
    double binSizeFactor = 1;
    binSizeFactor *= _mprime.getBinSize();
    binSizeFactor *= _thetaprime.getBinSize();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));
    }

    // This line runs once
    static std::array<fptype, 6> host_norms{{0, 0, 0, 0, 0, 0}};

    std::array<fptype, 6> current_host_norms{{_mprime.getLowerLimit(),
                                              _mprime.getUpperLimit(),
                                              static_cast<fptype>(_mprime.getNumBins()),
                                              _thetaprime.getLowerLimit(),
                                              _thetaprime.getUpperLimit(),
                                              static_cast<fptype>(_thetaprime.getNumBins())}};

    if(host_norms != current_host_norms) {
        host_norms = current_host_norms;
    }

    MEMCPY(dalitzNormRange, host_norms.data(), 6 * sizeof(fptype), cudaMemcpyHostToDevice);

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(decayInfo.resonances[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    forceRedoIntegrals = false;

    // Only do this bit if masses or widths have changed.
    thrust::constant_iterator<fptype *> arrayAddress(dalitzNormRange);
    thrust::counting_iterator<int> fakeEvtIndex(0);
    thrust::constant_iterator<fptype *> fakedataArray(dev_fake_event_array);
    thrust::constant_iterator<int> fakeeventSize(totalEventSize);

    for(int i = 0; i < n_res; ++i) {
        for(int j = 0; j < n_res ; ++j) {

            integrators_ff[i][j]->setDalitzIndex(getFunctionIndex());
            integrators_ff[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            integrators_ff[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());
            integrators_ff[i][j]->setNoEff();
            thrust::constant_iterator<int> effFunc(efficiencyFunction);
            fpcomplex dummy_ff(0, 0);
            thrust::plus<fpcomplex> complexSum_ff;
            (*(integrals_ff[i][j])) = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(fakeEvtIndex, fakedataArray, fakeeventSize, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(fakeEvtIndex + NumNormEvents, arrayAddress, fakeeventSize, effFunc)),
                *(integrators_ff[i][j]),
                dummy_ff,
                complexSum_ff);
            
        }

    }

    // End of time-consuming integrals.
    fpcomplex sumIntegral(0, 0);
    std::vector<std::vector<fptype>> AmpIntegral(n_res, std::vector<fptype>(n_res));

    for(unsigned int i = 0; i < n_res; ++i) {
        fpcomplex amplitude_i(host_parameters[parametersIdx + i * 2 + 1], host_parameters[parametersIdx + i * 2 + 2]);
        fpcomplex buffer(0., 0.);

        for(unsigned int j = 0; j < n_res; ++j) {
            fpcomplex amplitude_j(host_parameters[parametersIdx + j * 2 + 1],
                                  -host_parameters[parametersIdx + j * 2 + 2]);

            fptype norm_i = 1./((*(integrals[i][i])).real()/NumNormEvents);
            fptype norm_j = 1./((*(integrals[j][j])).real()/NumNormEvents);
            
            buffer = amplitude_i * amplitude_j * ((*(integrals[i][j]))/NumNormEvents)*sqrt(norm_i)*sqrt(norm_j);
    
            AmpIntegral[i][j] = buffer.real();
            sumIntegral += buffer;
        }
    }

    totalFF_integral = sumIntegral.real();

    for(int i = 0; i < n_res; i++) {
        for(int j = 0; j < n_res; j++) {
            AmpIntegral[i][j] /= totalFF_integral;
            AmpIntegral[i][j] *= 100;
        }
    }

    if(print) {
        std::cout << "Fit Fractions Matrix (%): \n";
        std::cout << "*Note: the order of diag FFs is equal to the order that which resonances are pushed into the "
                     "resonance vector. \n";
        Eigen::MatrixXd m(n_res, n_res);
        for(int i = 0; i < n_res; i++)
            m.row(i) = Eigen::Map<Eigen::VectorXd>(&AmpIntegral[i][0], n_res);

        std::cout << std::fixed << m << std::endl;
        fptype sumdiagffs = 0.;

        std::cout << "\n ";
        std::cout << "Diagonal Fit Fractions (%): \n";


        for(int i = 0; i < n_res; i++){
            auto name = decayInfo.resonances[i]->getName();
            std::cout  << name << "\t" << std::fixed << m(i, i) << '\n';
            sumdiagffs += m(i, i);
        }
        std::cout << "Sum of Diag FFs: " << sumdiagffs << "\n";
        std::cout << "\n";
    }

    return AmpIntegral;
}

__host__ auto Amp3BodySqDP::GenerateSig(unsigned int numEvents, int seed) -> std::
    tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h> {
    // Must configure our functions before any calculations!
    // setupObservables();
    // setIndices();

    initialize();

    // Defining phase space
    std::vector<mcbooster::GReal_t> masses{decayInfo.daug1Mass, decayInfo.daug2Mass, decayInfo.daug3Mass};
    mcbooster::PhaseSpace phsp(decayInfo.motherMass, masses, numEvents, generation_offset);

    if(seed != 0) {
        phsp.SetSeed(seed);
    } else {
        GOOFIT_INFO("Current generator seed {}, offset {}", phsp.GetSeed(), generation_offset);
    }

    // Generating numEvents events. Events are all generated inside the phase space with uniform distribution in
    // momentum space. Events must be weighted to have phase space distribution
    phsp.Generate(mcbooster::Vector4R(decayInfo.motherMass, 0.0, 0.0, 0.0));

    auto d1 = phsp.GetDaughters(0);
    auto d2 = phsp.GetDaughters(1);
    auto d3 = phsp.GetDaughters(2);

    mcbooster::ParticlesSet_d pset(3);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;

    auto SigGen_M12_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_M13_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_M23_d = mcbooster::RealVector_d(numEvents);

    mcbooster::VariableSet_d VarSet_d(3);
    VarSet_d[0] = &SigGen_M12_d;
    VarSet_d[1] = &SigGen_M23_d;
    VarSet_d[2] = &SigGen_M13_d;

    // Evaluating invariant masses for each event
    
    Dim2 eval = Dim2();
    mcbooster::EvaluateArray<Dim2>(eval, pset, VarSet_d);

    mcbooster::VariableSet_d GooVarSet_d(3);
    GooVarSet_d[0] = VarSet_d[0];
    GooVarSet_d[1] = VarSet_d[2];
    GooVarSet_d[2] = VarSet_d[1];

    auto h1 = new mcbooster::Particles_h(d1);
    auto h2 = new mcbooster::Particles_h(d2);
    auto h3 = new mcbooster::Particles_h(d3);

    mcbooster::ParticlesSet_h ParSet(3);
    ParSet[0] = h1;
    ParSet[1] = h2;
    ParSet[2] = h3;

    auto SigGen_M12_h = new mcbooster::RealVector_h(SigGen_M12_d);
    auto SigGen_M23_h = new mcbooster::RealVector_h(SigGen_M23_d);
    auto SigGen_M13_h = new mcbooster::RealVector_h(SigGen_M13_d);

    mcbooster::VariableSet_h VarSet(3);
    VarSet[0] = SigGen_M12_h;
    VarSet[1] = SigGen_M23_h;
    VarSet[2] = SigGen_M13_h;

    mcbooster::RealVector_d weights(phsp.GetWeights());
    phsp.FreeResources();

    auto DS = new mcbooster::RealVector_d(3 * numEvents);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 2; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 3);
        thrust::copy(GooVarSet_d[i]->begin(), GooVarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 2, DS->end(), 3);
    thrust::copy(eventNumber, eventNumber + numEvents, sr.begin());

    // Giving events to GooFit. Format of dev_evt_array must be (s12, s13, eventNumber). s23 is calculated automatically
    // in src/PDFs/physics/detail/SpecialResonanceCalculator.cu
    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(numEvents, 3);

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

    // evaluating amplitudes for generated events, amplitudes are incorporated in weights
    thrust::transform(
        results.begin(), results.end(), weights.begin(), weights.begin(), thrust::multiplies<mcbooster::GReal_t>());

    // Filing accept/reject flags for resonant distribution for each generated event
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
