#include "TddpThrustFunctor.hh"
#include <complex>
using std::complex; 

#include "TddpHelperFunctions.hh"
#include "TddpHelperFunctions.cxx"

const int resonanceSize = 8;   // Number of parameters to describe one resonance
const int resonanceOffset = 8; // Offset of the first resonance into the parameter index array 
// Offset is number of parameters, constant index, indices for tau, xmix, and ymix, index
// of resolution function, and finally number of resonances (not calculable from nP
// because we don't know what the efficiency and time resolution might need). Efficiency 
// and time-resolution parameters are after the resonance information. 

// The function of this array is to hold all the cached waves; specific 
// waves are recalculated when the corresponding resonance mass or width 
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone! 
__device__ WaveHolder* cWaves[10]; 

__device__ inline int parIndexFromResIndex (int resIndex) {
  return resonanceOffset + resIndex*resonanceSize; 
}

__device__ ThreeComplex device_Tddp_calcIntegrals (fptype m12, fptype m13, int res_i, int res_j, fptype* p, unsigned int* indices) {
  // For calculating Dalitz-plot integrals. What's needed is the products 
  // AiAj*, AiBj*, and BiBj*, where 
  // Ai = BW_i(x, y) + BW_i(y, x)
  // and Bi reverses the sign of the second BW. 
  // This function returns the above values at a single point. 
  // NB: Multiplication by efficiency is done by the calling function. 

  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3];  

  ThreeComplex ret; 
  //if (0 == threadIdx.x) cuPrintf("%f %f eval %i %f %f %f %f %f\n", m12, m13, indices[1], functorConstants[0], motherMass, daug1Mass, daug2Mass, daug3Mass);  
  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  //if (0 == threadIdx.x) cuPrintf("%f %f eval 2\n", m12, m13);  

  int parameter_i = parIndexFromResIndex(res_i);
  int parameter_j = parIndexFromResIndex(res_j);

  //fptype amp_real             = p[indices[parameter_i+0]];
  //fptype amp_imag             = p[indices[parameter_i+1]];
  fptype mass_i                 = p[indices[parameter_i+2]];
  fptype width_i                = p[indices[parameter_i+3]];
  unsigned int spin_i           = indices[parameter_i+4];
  unsigned int cyclic_index_i   = indices[parameter_i+5];
  unsigned int eval_type_i      = indices[parameter_i+6];

  // I don't use this anymore. We have wave A and wave B; 
  // it used to be the case that we calculated wave A as the sum of
  // flavour eigenstates plus CP-even states, and B as flavour plus CP-odd, 
  // because of the way the CP states cancel across the Dalitz plot. 
  // This was meant as an optimisation. However, since I need to calculate
  // the opposite-side numbers anyway to put them in the cache and 
  // use them for the normalisation integrals, it doesn't actually
  // save any time anymore. 
  //unsigned int resonance_type_i = indices[parameter_i+7];

  devcomplex<fptype> ai = getResonanceAmplitude(m12, m13, mass_i, width_i, spin_i, cyclic_index_i, eval_type_i, indices); 
  devcomplex<fptype> bi = getResonanceAmplitude(m13, m12, mass_i, width_i, spin_i, cyclic_index_i, eval_type_i, indices); 

  fptype mass_j                 = p[indices[parameter_j+2]];
  fptype width_j                = p[indices[parameter_j+3]];
  unsigned int spin_j           = indices[parameter_j+4];
  unsigned int cyclic_index_j   = indices[parameter_j+5];
  unsigned int eval_type_j      = indices[parameter_j+6];
  //unsigned int resonance_type_j = indices[parameter_j+7];
      
  devcomplex<fptype> aj = conj(getResonanceAmplitude(m12, m13, mass_j, width_j, spin_j, cyclic_index_j, eval_type_j, indices)); 
  devcomplex<fptype> bj = conj(getResonanceAmplitude(m13, m12, mass_j, width_j, spin_j, cyclic_index_j, eval_type_j, indices)); 

  //if ((3 == res_i) && (15 == res_j)) {
  //cuPrintf("Tddp integral 1: %f %f %f %f %i %i %i\n", m12, m13, mass_i, width_i, spin_i, cyclic_index_i, eval_type_i); 
  //cuPrintf("Tddp integral 2: %f %f %f %f %i %i %i\n", m12, m13, mass_j, width_j, spin_j, cyclic_index_j, eval_type_j); 
    //cuPrintf("Tddp integral 1: (%f, %f) (%f, %f) (%f, %f) (%f, %f) \n", m12, m13, (ai*aj).real, (ai*aj).imag, (ai*bj).real, (ai*bj).imag, (bi*bj).real, (bi*bj).imag); 
  //}
  
  ret = ThreeComplex((ai*aj).real, (ai*aj).imag, (ai*bj).real, (ai*bj).imag, (bi*bj).real, (bi*bj).imag);
  return ret; 
}

__device__ fptype device_Tddp (fptype* evt, fptype* p, unsigned int* indices) {
  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3]; 

  fptype m12 = evt[indices[4 + indices[0]]]; 
  fptype m13 = evt[indices[5 + indices[0]]];

  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return 0; 
  int evtNum = (int) FLOOR(0.5 + evt[indices[6 + indices[0]]]); 

  devcomplex<fptype> sumWavesA(0, 0);
  devcomplex<fptype> sumWavesB(0, 0); 
  devcomplex<fptype> sumRateAA(0, 0); 
  devcomplex<fptype> sumRateAB(0, 0); 
  devcomplex<fptype> sumRateBB(0, 0); 

  unsigned int numResonances = indices[6]; 
  unsigned int cacheToUse    = indices[7]; 

  for (int i = 0; i < numResonances; ++i) {
    int paramIndex  = parIndexFromResIndex(i);
    fptype amp_real = p[indices[paramIndex+0]];
    fptype amp_imag = p[indices[paramIndex+1]];

    devcomplex<fptype> matrixelement(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + i]),
				     thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + i])); 
    //if (evtNum == 43008) printf("Tddp: %i %i (%f %f) (%f %f) %i\n", evtNum, i, matrixelement.real, matrixelement.imag, amp_real, amp_imag, numResonances); 
    matrixelement.multiply(amp_real, amp_imag); 
    sumWavesA += matrixelement; 

#ifdef DEBUGSUMRATES
    if (25 > evtNum) {
      devcomplex<fptype> waveA_i(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + i]),
				 thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + i])); 
      devcomplex<fptype> waveB_i(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + i]),
				 thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + i])); 

      for (int j = 0; j < numResonances; ++j) {
	int paramIndex_j  = parIndexFromResIndex(j);
	fptype amp_real_j = p[indices[paramIndex_j+0]];
	fptype amp_imag_j = p[indices[paramIndex_j+1]];
	

	devcomplex<fptype> waveA_j(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + j]),
				   thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + j])); 

	devcomplex<fptype> waveB_j(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + j]),
				   thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + j])); 
	devcomplex<fptype> amps(amp_real, -amp_imag);
	amps.multiply(amp_real_j, amp_imag_j); 

	devcomplex<fptype> rateAA = conj(waveA_i)*waveA_j*amps;
	devcomplex<fptype> rateAB = conj(waveA_i)*waveB_j*amps;
	devcomplex<fptype> rateBB = conj(waveB_i)*waveB_j*amps;

	//printf("Wave (%i, %i): (%f, %f) (%f, %f) (%f, %f) \n", i, j, rateAA.real, rateAA.imag, rateAB.real, rateAB.imag, rateBB.real, rateBB.imag); 
	
	sumRateAA += rateAA;
	sumRateAB += rateAB;
	sumRateBB += rateBB;
      }
      waveA_i.multiply(amp_real, amp_imag);
      waveB_i.multiply(amp_real, amp_imag);
      //printf("WaveA,B %i (%f, %f), (%f, %f)\n", i, waveA_i.real, waveA_i.imag, waveB_i.real, waveB_i.imag); 
    }
#endif
    

    matrixelement = devcomplex<fptype>(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + i]),
				       thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + i])); 
    matrixelement.multiply(amp_real, amp_imag); 
    sumWavesB += matrixelement; 
  } 

  fptype _tau     = p[indices[2]];
  fptype _xmixing = p[indices[3]];
  fptype _ymixing = p[indices[4]];
  
  fptype _time    = evt[indices[2 + indices[0]]];
  fptype _sigma   = evt[indices[3 + indices[0]]];

  //if ((gpuDebug & 1) && (0 == blockIdx.x) && (0 == threadIdx.x)) 
  //if (0 == evtNum) printf("TDDP: (%f, %f) (%f, %f)\n", sumWavesA.real, sumWavesA.imag, sumWavesB.real, sumWavesB.imag);
  //printf("TDDP: %f %f %f %f | %f %f %i\n", m12, m13, _time, _sigma, _xmixing, _tau, evtNum); 


  /*
  fptype ret = 0; 
  ret += (norm2(sumWavesA) + norm2(sumWavesB))*COSH(_ymixing * _time);
  ret += (norm2(sumWavesA) - norm2(sumWavesB))*COS (_xmixing * _time);
  sumWavesA *= conj(sumWavesB); 
  ret -= 2*sumWavesA.real * SINH(_ymixing * _time);
  ret -= 2*sumWavesA.imag * SIN (_xmixing * _time); // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B. 
  ret *= EXP(-_time); 
  */

  fptype term1 = norm2(sumWavesA) + norm2(sumWavesB);
  fptype term2 = norm2(sumWavesA) - norm2(sumWavesB);
  sumWavesA *= conj(sumWavesB); 
  //printf("(%i, %i) TDDP: %f %f %f %f %f %f %f\n", blockIdx.x, threadIdx.x, term1, term2, sumWavesA.real, sumWavesA.imag, m12, m13, _tau);

  // Cannot use callFunction on resolution function. 
  int effFunctionIdx = parIndexFromResIndex(numResonances); 
  fptype ret = (*(reinterpret_cast<device_resfunction_ptr>(device_function_table[indices[5]])))(term1, term2, sumWavesA.real, sumWavesA.imag,
												_tau, _time, _xmixing, _ymixing, _sigma, 
												p, &(indices[2 + effFunctionIdx])); 
  
  //if ((evtNum < 1) && (gpuDebug & 1)) internalDebug = 2;
  //fptype eff = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[effFunctionIdx]])))(evt, p, paramIndices + indices[effFunctionIdx + 1]);
  fptype eff = callFunction(evt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  //internalDebug = 0; 
  ret *= eff;

  //if ((gpuDebug & 1) && (4000 > evtNum)) {
    //printf("FULLPRINT: %i: %f %f\n", evtNum, ret*normalisationFactors[(indices - paramIndices)], eff);
    //printf("FULLPRINT: %i: %f %f (%f %f %f %f)\n", evtNum, ret, eff, m12, m13, _time, _sigma); 
  //} 



  //if ((evtNum < 50) && (gpuDebug & 1)) {
  //if ((gpuDebug & 1) && (0 == threadIdx.x)) {
  //if ((gpuDebug & 1) && (180 == evtNum)) {
  //if ((0 == threadIdx.x) && (0 == blockIdx.x) && (gpuDebug & 1)) {
  //sumRateAA *= eff;
  //sumRateAB *= eff;
  //sumRateBB *= eff;
  //printf("Signal1 %i (%f, %f) (%f, %f) (%f, %f)\n", evtNum, sumRateAA.real, sumRateAA.imag, sumRateAB.real, sumRateAB.imag, sumRateBB.real, sumRateBB.imag); 
  //printf("Signal1 %i (%f, %f) (%f, %f) (%f, %f %f)\n", evtNum, term1, term2, sumWavesA.real, sumWavesA.imag, _tau, _xmixing, _ymixing);
  //fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 
  //printf("Signal2 %i (%f, %f, %f) %f, %f | %f %f %f\n", evtNum, m12, m13, m23, _time, _sigma, eff, ret, normalisationFactors[(indices - paramIndices)]);
  //printf("Signal3 %f %f %f %f %f %f %f %f\n", cudaArray[indices[effFunctionIdx+1]+1], cudaArray[indices[effFunctionIdx+1]+2], cudaArray[indices[effFunctionIdx+1]+3], 
  //cudaArray[indices[effFunctionIdx+1]+4], cudaArray[indices[effFunctionIdx+1]+5], cudaArray[indices[effFunctionIdx+1]+6], 
  //cudaArray[indices[effFunctionIdx+1]+7], cudaArray[indices[effFunctionIdx+1]+8]); 
  //}

  //printf("(%i, %i) TDDP: %f %f %f %f %f %i %f\n", blockIdx.x, threadIdx.x, _time, _sigma, m12, m13, term1, evtNum, ret);
  //if ((gpuDebug & 1) && (isnan(ret)))
  //printf("(%i, %i) TDDP: %f %f %f %f %i %i %f\n", blockIdx.x, threadIdx.x, _time, _sigma, m12, m13, evtNum, indices[6 + indices[0]], evt[indices[6 + indices[0]]]);
  //if ((gpuDebug & 1) && (isnan(ret)))
  //printf("(%i, %i) TDDP: %f %f %f %f %f %f %f\n", blockIdx.x, threadIdx.x, term1, term2, sumWavesA.real, sumWavesA.imag, _xmixing, _ymixing, _tau);

  return ret; 
}

__device__ device_function_ptr ptr_to_Tddp = device_Tddp; 

__host__ TddpThrustFunctor::TddpThrustFunctor (std::string n, Variable* _dtime, Variable* _sigmat, Variable* m12, Variable* m13, Variable* eventNumber, DecayInfo* decay, MixingTimeResolution* r, ThrustPdfFunctor* efficiency) 
  : ThrustPdfFunctor(_dtime, n) 
  , decayInfo(decay)
  , _m12(m12)
  , _m13(m13)
  , dalitzNormRange(0)
  , cachedWaves(0) 
  , integrals(0)
  , resolution(r) 
  , forceRedoIntegrals(true)
  , totalEventSize(5) // Default 5 = m12, m13, time, sigma_t, evtNum 
  , cacheToUse(0) 
  , integrators(0)
  , calculators(0) 
{
  // NB, _dtime already registered!
  registerObservable(_sigmat);
  registerObservable(_m12);
  registerObservable(_m13);
  registerObservable(eventNumber); 

  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(5)); 
  fptype decayConstants[5];
  decayConstants[0] = decayInfo->motherMass;
  decayConstants[1] = decayInfo->daug1Mass;
  decayConstants[2] = decayInfo->daug2Mass;
  decayConstants[3] = decayInfo->daug3Mass;
  decayConstants[4] = decayInfo->meson_radius;
  cudaMemcpyToSymbol(functorConstants, decayConstants, 5*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);  
  
  pindices.push_back(registerParameter(decayInfo->_tau));
  pindices.push_back(registerParameter(decayInfo->_xmixing));
  pindices.push_back(registerParameter(decayInfo->_ymixing));
  assert(resolution->getDeviceFunction() >= 0); 
  pindices.push_back((unsigned int) resolution->getDeviceFunction()); 
  pindices.push_back(decayInfo->resonances.size()); 

  static int cacheCount = 0; 
  cacheToUse = cacheCount++; 
  pindices.push_back(cacheToUse); 

  for (std::vector<ResonanceInfo*>::iterator res = decayInfo->resonances.begin(); res != decayInfo->resonances.end(); ++res) {
    pindices.push_back(registerParameter((*res)->amp_real));
    pindices.push_back(registerParameter((*res)->amp_imag));
    pindices.push_back(registerParameter((*res)->mass));
    pindices.push_back(registerParameter((*res)->width));
    pindices.push_back((*res)->spin);
    pindices.push_back((*res)->cyclic_index);
    pindices.push_back((*res)->eval_type);
    pindices.push_back((*res)->resonance_type); 
  }

  pindices.push_back(efficiency->getFunctionIndex());
  pindices.push_back(efficiency->getParameterIndex());
  components.push_back(efficiency); 

  resolution->createParameters(pindices, this); 
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_Tddp, sizeof(void*));
  initialise(pindices);

  redoIntegral = new bool[decayInfo->resonances.size()];
  cachedMasses = new fptype[decayInfo->resonances.size()];
  cachedWidths = new fptype[decayInfo->resonances.size()];
  integrals    = new ThreeComplex**[decayInfo->resonances.size()];
  integrators  = new SpecialDalitzIntegrator**[decayInfo->resonances.size()];
  calculators  = new SpecialWaveCalculator*[decayInfo->resonances.size()];

  for (int i = 0; i < decayInfo->resonances.size(); ++i) {
    redoIntegral[i] = true;
    cachedMasses[i] = -1;
    cachedWidths[i] = -1; 
    integrators[i]  = new SpecialDalitzIntegrator*[decayInfo->resonances.size()];
    calculators[i]  = new SpecialWaveCalculator(parameters, i); 
    integrals[i]    = new ThreeComplex*[decayInfo->resonances.size()];
    
    for (int j = 0; j < decayInfo->resonances.size(); ++j) {
      integrals[i][j]   = new ThreeComplex(0, 0, 0, 0, 0, 0);
      integrators[i][j] = new SpecialDalitzIntegrator(parameters, i, j); 
    }
  }

  addSpecialMask(FunctorBase::ForceSeparateNorm); 
}

__host__ void TddpThrustFunctor::setDataSize (unsigned int dataSize, unsigned int evtSize) {
  // Default 5 is m12, m13, time, sigma_t, evtNum
  totalEventSize = evtSize;
  assert(totalEventSize >= 5); 

  if (cachedWaves) {
    delete cachedWaves;
  }

  numEntries = dataSize; 
  cachedWaves = new thrust::device_vector<WaveHolder>(dataSize*decayInfo->resonances.size());
  void* dummy = thrust::raw_pointer_cast(cachedWaves->data()); 
  cudaMemcpyToSymbol(cWaves, &dummy, sizeof(WaveHolder*), cacheToUse*sizeof(WaveHolder*)); 
  setForceIntegrals(); 
}

__host__ fptype TddpThrustFunctor::normalise () const {
  recursiveSetNormalisation(1); // Not going to normalise efficiency, 
  // so set normalisation factor to 1 so it doesn't get multiplied by zero. 
  // Copy at this time to ensure that the SpecialWaveCalculators, which need the efficiency, 
  // don't get zeroes through multiplying by the normFactor. 
  cudaMemcpyToSymbol(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  //std::cout << "TDDP normalisation " << getName() << std::endl;

  int totalBins = _m12->numbins * _m13->numbins;
  if (!dalitzNormRange) {
    cudaMalloc((void**) &dalitzNormRange, 6*sizeof(fptype));
  
    fptype* host_norms = new fptype[6];
    host_norms[0] = _m12->lowerlimit;
    host_norms[1] = _m12->upperlimit;
    host_norms[2] = _m12->numbins;
    host_norms[3] = _m13->lowerlimit;
    host_norms[4] = _m13->upperlimit;
    host_norms[5] = _m13->numbins;
    cudaMemcpy(dalitzNormRange, host_norms, 6*sizeof(fptype), cudaMemcpyHostToDevice);
    delete[] host_norms; 
  }

  for (unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
    int param_i = parameters + resonanceOffset + resonanceSize*i; 
    redoIntegral[i] = forceRedoIntegrals; 
    if ((host_params[host_indices[param_i + 2]] == cachedMasses[i]) && (host_params[host_indices[param_i + 3]] == cachedWidths[i])) continue;
    redoIntegral[i] = true; 
    cachedMasses[i] = host_params[host_indices[param_i + 2]];
    cachedWidths[i] = host_params[host_indices[param_i + 3]];
  }
  forceRedoIntegrals = false; 

  // Only do this bit if masses or widths have changed.  
  thrust::constant_iterator<fptype*> arrayAddress(dalitzNormRange); 
  thrust::counting_iterator<int> binIndex(0); 

  // NB, SpecialWaveCalculator assumes that fit is unbinned! 
  // And it needs to know the total event size, not just observables
  // for this particular PDF component. 
  thrust::constant_iterator<fptype*> dataArray(cudaDataArray); 
  thrust::constant_iterator<int> eventSize(totalEventSize);
  thrust::counting_iterator<int> eventIndex(0); 

  static int normCall = 0; 
  normCall++; 

  for (int i = 0; i < decayInfo->resonances.size(); ++i) {
    if (redoIntegral[i]) {
      
      thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
			thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
			strided_range<thrust::device_vector<WaveHolder>::iterator>(cachedWaves->begin() + i, cachedWaves->end(), decayInfo->resonances.size()).begin(), 
			*(calculators[i]));
      //std::cout << "Integral for resonance " << i << " " << numEntries << " " << totalEventSize << std::endl; 
      /*
      WaveHolder hostWH = (*cachedWaves)[i]; 
      
      std::cout << "Resonance " << i << " "
		<< normCall << " : (" 
		<< hostWH.waveAreal << ", " << hostWH.waveAimag << ") ("
		<< hostWH.waveBreal << ", " << hostWH.waveBimag << ")\n";
      */
    }
    
    // Possibly this can be done more efficiently by exploiting symmetry? 
    for (int j = 0; j < decayInfo->resonances.size(); ++j) {
      if ((!redoIntegral[i]) && (!redoIntegral[j])) continue; 
      ThreeComplex dummy(0, 0, 0, 0, 0, 0);
      SpecialComplexSum complexSum; 
      (*(integrals[i][j])) = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
						      thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
						      *(integrators[i][j]), 
						      dummy, 
						      complexSum); 
    }
  }      

  // End of time-consuming integrals. 

  complex<fptype> integralA_2(0, 0);
  complex<fptype> integralB_2(0, 0);
  complex<fptype> integralABs(0, 0);
  for (unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
    int param_i = parameters + resonanceOffset + resonanceSize*i; 
    complex<fptype> amplitude_i(host_params[host_indices[param_i]], host_params[host_indices[param_i + 1]]);

    for (unsigned int j = 0; j < decayInfo->resonances.size(); ++j) {
      int param_j = parameters + resonanceOffset + resonanceSize*j; 
      complex<fptype> amplitude_j(host_params[host_indices[param_j]], -host_params[host_indices[param_j + 1]]); // Notice complex conjugation

      integralA_2 += (amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j])))); 
      integralABs += (amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j]))));
      integralB_2 += (amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j]))));

      /*
      if (cpuDebug & 1) {
	int idx = i * decayInfo->resonances.size() + j;     
	if (0 == host_callnumber) std::cout << "Integral contribution " << i << ", " << j << " " << idx << " : "
					    << amplitude_i << " "
					    << amplitude_j << " (" 
					    << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j])))) << ", "
					    << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j])))) << ") ("
					    << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j])))) << ", "
					    << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j])))) << ") ("
					    << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j])))) << ", "
					    << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j])))) << ") "
					    << thrust::get<0>(*(integrals[i][j])) << ", "
					    << thrust::get<1>(*(integrals[i][j])) << ") ("
					    << thrust::get<2>(*(integrals[i][j])) << ", "
					    << thrust::get<3>(*(integrals[i][j])) << ") ("
					    << thrust::get<4>(*(integrals[i][j])) << ", "
					    << thrust::get<5>(*(integrals[i][j])) << ") ("
					    << real(integralA_2) << ", " << imag(integralA_2) << ") "
					    << std::endl; 
      }
      */

    }
  }

  double dalitzIntegralOne = real(integralA_2); // Notice that this is already the abs2, so it's real by construction; but the compiler doesn't know that. 
  double dalitzIntegralTwo = real(integralB_2);
  double dalitzIntegralThr = real(integralABs);
  double dalitzIntegralFou = imag(integralABs); 
  
  fptype tau     = host_params[host_indices[parameters + 2]];
  fptype xmixing = host_params[host_indices[parameters + 3]];
  fptype ymixing = host_params[host_indices[parameters + 4]];

  fptype ret = resolution->normalisation(dalitzIntegralOne, dalitzIntegralTwo, dalitzIntegralThr, dalitzIntegralFou, tau, xmixing, ymixing); 

  /*
  fptype timeIntegralOne = tau / (1 - ymixing*ymixing); 
  fptype timeIntegralTwo = tau / (1 + xmixing*xmixing);
  fptype timeIntegralThr = ymixing * timeIntegralOne;
  fptype timeIntegralFou = xmixing * timeIntegralTwo;
       
  double ret = timeIntegralOne * (dalitzIntegralOne + dalitzIntegralTwo);
  ret       += timeIntegralTwo * (dalitzIntegralOne - dalitzIntegralTwo);
  ret       -= 2*timeIntegralThr * dalitzIntegralThr;
  ret       -= 2*timeIntegralFou * dalitzIntegralFou;
  */
  double binSizeFactor = 1;
  binSizeFactor *= ((_m12->upperlimit - _m12->lowerlimit) / _m12->numbins);
  binSizeFactor *= ((_m13->upperlimit - _m13->lowerlimit) / _m13->numbins);
  ret *= binSizeFactor;

  //cudaPrintfDisplay(stdout, true);
  //cudaPrintfEnd();

#ifdef CUDAPRINT 
  std::cout << "Tddp dalitz integrals: " 
	    << dalitzIntegralOne * binSizeFactor << " "
	    << dalitzIntegralTwo * binSizeFactor << " "
	    << dalitzIntegralThr * binSizeFactor << " "
	    << dalitzIntegralFou * binSizeFactor << " | "
	    << binSizeFactor << " " 
	    << tau << " " << xmixing << " " << ymixing << " "
	    << ret << " " 
	    << std::endl; 
#endif
    
  host_normalisation[parameters] = 1.0/ret;
  //std::cout << "End of TDDP normalisation: " << ret << " " << host_normalisation[parameters] << " " << binSizeFactor << std::endl; 
  return (fptype) ret; 
}
//#endif 

SpecialDalitzIntegrator::SpecialDalitzIntegrator (int pIdx, unsigned int ri, unsigned int rj) 
  : resonance_i(ri)
  , resonance_j(rj)
  , parameters(pIdx) 
{}

__device__ ThreeComplex SpecialDalitzIntegrator::operator () (thrust::tuple<int, fptype*> t) const {
  // Bin index, base address [lower, upper, numbins] 
  // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
  // that event size is two, and that the function to call is dev_Tddp_calcIntegrals.

  int globalBinNumber  = thrust::get<0>(t);
  fptype lowerBoundM12 = thrust::get<1>(t)[0];
  fptype upperBoundM12 = thrust::get<1>(t)[1];  
  int numBinsM12       = (int) FLOOR(thrust::get<1>(t)[2] + 0.5); 
  int binNumberM12     = globalBinNumber % numBinsM12;
  fptype binCenterM12  = upperBoundM12 - lowerBoundM12;
  binCenterM12        /= numBinsM12;
  binCenterM12        *= (binNumberM12 + 0.5); 
  binCenterM12        += lowerBoundM12; 

  globalBinNumber     /= numBinsM12; 
  fptype lowerBoundM13 = thrust::get<1>(t)[3];
  fptype upperBoundM13 = thrust::get<1>(t)[4];  
  int numBinsM13       = (int) FLOOR(thrust::get<1>(t)[5] + 0.5); 
  fptype binCenterM13  = upperBoundM13 - lowerBoundM13;
  binCenterM13        /= numBinsM13;
  binCenterM13        *= (globalBinNumber + 0.5); 
  binCenterM13        += lowerBoundM13; 

  //if (0 == threadIdx.x) cuPrintf("%i %i %i %f %f operator\n", thrust::get<0>(t), thrust::get<0>(t) % numBinsM12, globalBinNumber, binCenterM12, binCenterM13);
  unsigned int* indices = paramIndices + parameters;   
  ThreeComplex ret = device_Tddp_calcIntegrals(binCenterM12, binCenterM13, resonance_i, resonance_j, cudaArray, indices); 

  fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an event-weighted fit. 
  fakeEvt[indices[indices[0] + 2 + 2]] = binCenterM12;
  fakeEvt[indices[indices[0] + 2 + 3]] = binCenterM13;
  unsigned int numResonances = indices[6]; 
  int effFunctionIdx = parIndexFromResIndex(numResonances); 
  //if (thrust::get<0>(t) == 19840) {internalDebug1 = blockIdx.x; internalDebug2 = threadIdx.x;}
  //fptype eff = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[effFunctionIdx]])))(fakeEvt, cudaArray, paramIndices + indices[effFunctionIdx + 1]);
  fptype eff = callFunction(fakeEvt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  //if (thrust::get<0>(t) == 19840) {
  //internalDebug1 = -1; 
  //internalDebug2 = -1;
    //printf("Efficiency: %i %f %f %f %i\n", thrust::get<0>(t), binCenterM12, binCenterM13, eff, effFunctionIdx);
    //printf("Efficiency: %f %f %f %f %f %i %i\n", fakeEvt[0], fakeEvt[1], fakeEvt[2], fakeEvt[3], fakeEvt[4], indices[indices[0] + 2 + 2], indices[indices[0] + 2 + 3]); 
  //}

  // Multiplication by eff, not sqrt(eff), is correct:
  // These complex numbers will not be squared when they
  // go into the integrals. They've been squared already,
  // as it were. 
  thrust::get<0>(ret) *= eff;
  thrust::get<1>(ret) *= eff;
  thrust::get<2>(ret) *= eff;
  thrust::get<3>(ret) *= eff;
  thrust::get<4>(ret) *= eff;
  thrust::get<5>(ret) *= eff;
  return ret; 
}

SpecialWaveCalculator::SpecialWaveCalculator (int pIdx, unsigned int res_idx) 
  : resonance_i(res_idx)
  , parameters(pIdx)
{}

__device__ WaveHolder SpecialWaveCalculator::operator () (thrust::tuple<int, fptype*, int> t) const {
  // Calculates the BW values for a specific resonance. 

  WaveHolder ret;

  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t)); 

  unsigned int* indices = paramIndices + parameters;   // Jump to TDDP position within parameters array
  fptype m12 = evt[indices[4 + indices[0]]]; 
  fptype m13 = evt[indices[5 + indices[0]]];

  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3];  

  //printf("(%i, %i) Special: %i %i %f %f\n", blockIdx.x, threadIdx.x, evtNum, thrust::get<2>(t), m12, m13); 
  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;

  int parameter_i = parIndexFromResIndex(resonance_i); // Find position of this resonance relative to TDDP start 

  // Now extract actual parameters using indices found above. Notice double indirection: 'parameters'  into paramIndices, then paramIndices content into cudaArray. 
  fptype mass_i                 = cudaArray[indices[parameter_i+2]];
  fptype width_i                = cudaArray[indices[parameter_i+3]];
  unsigned int spin_i           = indices[parameter_i+4]; // Bit of a kludge to allow constant integer `parameters' without casting from functorConstants. 
  unsigned int cyclic_index_i   = indices[parameter_i+5];
  unsigned int eval_type_i      = indices[parameter_i+6];

  devcomplex<fptype> ai = getResonanceAmplitude(m12, m13, mass_i, width_i, spin_i, cyclic_index_i, eval_type_i, indices); 
  devcomplex<fptype> bi = getResonanceAmplitude(m13, m12, mass_i, width_i, spin_i, cyclic_index_i, eval_type_i, indices); 

  //if (evtNum < 1) cuPrintf("Wave %i %i: (%f, %f) (%f, %f)\n", evtNum, resonance_i, ai.real, ai.imag, bi.real, bi.imag); 

  thrust::get<0>(ret) = ai.real;
  thrust::get<1>(ret) = ai.imag;
  thrust::get<2>(ret) = bi.real;
  thrust::get<3>(ret) = bi.imag; 

  return ret; 
}

