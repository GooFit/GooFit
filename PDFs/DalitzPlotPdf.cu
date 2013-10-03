#include "DalitzPlotPdf.hh"
#include <complex>
using std::complex; 

const int resonanceOffset_DP = 4; // Offset of the first resonance into the parameter index array 
// Offset is number of parameters, constant index, number of resonances (not calculable 
// from nP because we don't know what the efficiency might need), and cache index. Efficiency 
// parameters are after the resonance information. 

// The function of this array is to hold all the cached waves; specific 
// waves are recalculated when the corresponding resonance mass or width 
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone! 
MEM_DEVICE devcomplex<fptype>* cResonances[10]; 

EXEC_TARGET inline int parIndexFromResIndex_DP (int resIndex) {
  return resonanceOffset_DP + resIndex*resonanceSize; 
}

EXEC_TARGET devcomplex<fptype> device_DalitzPlot_calcIntegrals (fptype m12, fptype m13, int res_i, int res_j, fptype* p, unsigned int* indices) {
  // Calculates BW_i(m12, m13) * BW_j^*(m12, m13). 
  // This calculation is in a separate function so
  // it can be cached. Note that this function expects
  // to be called on a normalisation grid, not on 
  // observed points, that's why it doesn't use 
  // cResonances. No need to cache the values at individual
  // grid points - we only care about totals. 
  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3];  

  devcomplex<fptype> ret; 
  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 

  int parameter_i = parIndexFromResIndex_DP(res_i);
  unsigned int functn_i = indices[parameter_i+2];
  unsigned int params_i = indices[parameter_i+3];
  ret = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);

  int parameter_j = parIndexFromResIndex_DP(res_j);
  unsigned int functn_j = indices[parameter_j+2];
  unsigned int params_j = indices[parameter_j+3];
  ret *= conj(getResonanceAmplitude(m12, m13, m23, functn_j, params_j));

  return ret; 
}

EXEC_TARGET fptype device_DalitzPlot (fptype* evt, fptype* p, unsigned int* indices) {
  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3]; 

  fptype m12 = evt[indices[2 + indices[0]]]; 
  fptype m13 = evt[indices[3 + indices[0]]];

  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return 0; 
  int evtNum = (int) FLOOR(0.5 + evt[indices[4 + indices[0]]]); 

  devcomplex<fptype> totalAmp(0, 0);
  unsigned int numResonances = indices[2]; 
  unsigned int cacheToUse    = indices[3]; 

  for (int i = 0; i < numResonances; ++i) {
    int paramIndex  = parIndexFromResIndex_DP(i);
    fptype amp_real = p[indices[paramIndex+0]];
    fptype amp_imag = p[indices[paramIndex+1]];

    devcomplex<fptype> matrixelement((cResonances[cacheToUse][evtNum*numResonances + i]).real,
				     (cResonances[cacheToUse][evtNum*numResonances + i]).imag); 
    matrixelement.multiply(amp_real, amp_imag); 
    totalAmp += matrixelement; 
  } 
   
  fptype ret = norm2(totalAmp); 
  int effFunctionIdx = parIndexFromResIndex_DP(numResonances); 
  fptype eff = callFunction(evt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  ret *= eff;

  //printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real, totalAmp.imag); 

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_DalitzPlot = device_DalitzPlot; 

__host__ DalitzPlotPdf::DalitzPlotPdf (std::string n, 
							   Variable* m12, 
							   Variable* m13, 
							   Variable* eventNumber, 
							   DecayInfo* decay, 
							   GooPdf* efficiency)
  : GooPdf(0, n) 
  , decayInfo(decay)
  , _m12(m12)
  , _m13(m13)
  , dalitzNormRange(0)
  , cachedWaves(0) 
  , integrals(0)
  , forceRedoIntegrals(true)
  , totalEventSize(3) // Default 3 = m12, m13, evtNum 
  , cacheToUse(0) 
  , integrators(0)
  , calculators(0) 
{
  registerObservable(_m12);
  registerObservable(_m13);
  registerObservable(eventNumber); 

  fptype decayConstants[5];
  
  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(5)); 
  decayConstants[0] = decayInfo->motherMass;
  decayConstants[1] = decayInfo->daug1Mass;
  decayConstants[2] = decayInfo->daug2Mass;
  decayConstants[3] = decayInfo->daug3Mass;
  decayConstants[4] = decayInfo->meson_radius;
  MEMCPY_TO_SYMBOL(functorConstants, decayConstants, 5*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);  

  pindices.push_back(decayInfo->resonances.size()); 
  static int cacheCount = 0; 
  cacheToUse = cacheCount++; 
  pindices.push_back(cacheToUse); 

  for (std::vector<ResonancePdf*>::iterator res = decayInfo->resonances.begin(); res != decayInfo->resonances.end(); ++res) {
    pindices.push_back(registerParameter((*res)->amp_real));
    pindices.push_back(registerParameter((*res)->amp_imag));
    pindices.push_back((*res)->getFunctionIndex());
    pindices.push_back((*res)->getParameterIndex());
    (*res)->setConstantIndex(cIndex); 
    components.push_back(*res);
  }

  pindices.push_back(efficiency->getFunctionIndex());
  pindices.push_back(efficiency->getParameterIndex());
  components.push_back(efficiency); 

  GET_FUNCTION_ADDR(ptr_to_DalitzPlot);
  initialise(pindices);

  redoIntegral = new bool[decayInfo->resonances.size()];
  cachedMasses = new fptype[decayInfo->resonances.size()];
  cachedWidths = new fptype[decayInfo->resonances.size()];
  integrals    = new devcomplex<fptype>**[decayInfo->resonances.size()];
  integrators  = new SpecialResonanceIntegrator**[decayInfo->resonances.size()];
  calculators  = new SpecialResonanceCalculator*[decayInfo->resonances.size()];

  for (int i = 0; i < decayInfo->resonances.size(); ++i) {
    redoIntegral[i] = true;
    cachedMasses[i] = -1;
    cachedWidths[i] = -1; 
    integrators[i]  = new SpecialResonanceIntegrator*[decayInfo->resonances.size()];
    calculators[i]  = new SpecialResonanceCalculator(parameters, i); 
    integrals[i]    = new devcomplex<fptype>*[decayInfo->resonances.size()];
    
    for (int j = 0; j < decayInfo->resonances.size(); ++j) {
      integrals[i][j]   = new devcomplex<fptype>(0, 0); 
      integrators[i][j] = new SpecialResonanceIntegrator(parameters, i, j); 
    }
  }

  addSpecialMask(PdfBase::ForceSeparateNorm); 
}

__host__ void DalitzPlotPdf::setDataSize (unsigned int dataSize, unsigned int evtSize) {
  // Default 3 is m12, m13, evtNum
  totalEventSize = evtSize;
  assert(totalEventSize >= 3); 

  if (cachedWaves) delete cachedWaves;

  numEntries = dataSize; 
  cachedWaves = new DEVICE_VECTOR<devcomplex<fptype> >(dataSize*decayInfo->resonances.size());
  void* dummy = thrust::raw_pointer_cast(cachedWaves->data()); 
  MEMCPY_TO_SYMBOL(cResonances, &dummy, sizeof(devcomplex<fptype>*), cacheToUse*sizeof(devcomplex<fptype>*), cudaMemcpyHostToDevice); 
  setForceIntegrals(); 
}

__host__ fptype DalitzPlotPdf::normalise () const {
  recursiveSetNormalisation(1); // Not going to normalise efficiency, 
  // so set normalisation factor to 1 so it doesn't get multiplied by zero. 
  // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency, 
  // don't get zeroes through multiplying by the normFactor. 
  MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice); 

  int totalBins = _m12->numbins * _m13->numbins;
  if (!dalitzNormRange) {
    gooMalloc((void**) &dalitzNormRange, 6*sizeof(fptype));
  
    fptype* host_norms = new fptype[6];
    host_norms[0] = _m12->lowerlimit;
    host_norms[1] = _m12->upperlimit;
    host_norms[2] = _m12->numbins;
    host_norms[3] = _m13->lowerlimit;
    host_norms[4] = _m13->upperlimit;
    host_norms[5] = _m13->numbins;
    MEMCPY(dalitzNormRange, host_norms, 6*sizeof(fptype), cudaMemcpyHostToDevice);
    delete[] host_norms; 
  }

  for (unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
    redoIntegral[i] = forceRedoIntegrals; 
    if (!(decayInfo->resonances[i]->parametersChanged())) continue;
    redoIntegral[i] = true; 
    decayInfo->resonances[i]->storeParameters();
  }
  forceRedoIntegrals = false; 

  // Only do this bit if masses or widths have changed.  
  thrust::constant_iterator<fptype*> arrayAddress(dalitzNormRange); 
  thrust::counting_iterator<int> binIndex(0); 

  // NB, SpecialResonanceCalculator assumes that fit is unbinned! 
  // And it needs to know the total event size, not just observables
  // for this particular PDF component. 
  thrust::constant_iterator<fptype*> dataArray(dev_event_array); 
  thrust::constant_iterator<int> eventSize(totalEventSize);
  thrust::counting_iterator<int> eventIndex(0); 

  for (int i = 0; i < decayInfo->resonances.size(); ++i) {
    if (redoIntegral[i]) {
      thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
			thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
			strided_range<DEVICE_VECTOR<devcomplex<fptype> >::iterator>(cachedWaves->begin() + i, 
										    cachedWaves->end(), 
										    decayInfo->resonances.size()).begin(), 
			*(calculators[i]));
    }
    
    // Possibly this can be done more efficiently by exploiting symmetry? 
    for (int j = 0; j < decayInfo->resonances.size(); ++j) {
      if ((!redoIntegral[i]) && (!redoIntegral[j])) continue; 
      devcomplex<fptype> dummy(0, 0);
      thrust::plus<devcomplex<fptype> > complexSum; 
      (*(integrals[i][j])) = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
						      thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
						      *(integrators[i][j]), 
						      dummy, 
						      complexSum); 
    }
  }      

  // End of time-consuming integrals. 
  complex<fptype> sumIntegral(0, 0);
  for (unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
    int param_i = parameters + resonanceOffset_DP + resonanceSize*i; 
    complex<fptype> amplitude_i(host_params[host_indices[param_i]], host_params[host_indices[param_i + 1]]);
    for (unsigned int j = 0; j < decayInfo->resonances.size(); ++j) {
      int param_j = parameters + resonanceOffset_DP + resonanceSize*j; 
      complex<fptype> amplitude_j(host_params[host_indices[param_j]], -host_params[host_indices[param_j + 1]]); 
      // Notice complex conjugation

      sumIntegral += (amplitude_i * amplitude_j * complex<fptype>((*(integrals[i][j])).real, (*(integrals[i][j])).imag)); 
    }
  }

  fptype ret = real(sumIntegral); // That complex number is a square, so it's fully real
  double binSizeFactor = 1;
  binSizeFactor *= ((_m12->upperlimit - _m12->lowerlimit) / _m12->numbins);
  binSizeFactor *= ((_m13->upperlimit - _m13->lowerlimit) / _m13->numbins);
  ret *= binSizeFactor;

  host_normalisation[parameters] = 1.0/ret;
  return (fptype) ret; 
}

SpecialResonanceIntegrator::SpecialResonanceIntegrator (int pIdx, unsigned int ri, unsigned int rj) 
  : resonance_i(ri)
  , resonance_j(rj)
  , parameters(pIdx) 
{}

EXEC_TARGET devcomplex<fptype> SpecialResonanceIntegrator::operator () (thrust::tuple<int, fptype*> t) const {
  // Bin index, base address [lower, upper, numbins] 
  // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
  // that event size is two, and that the function to call is dev_DalitzPlot_calcIntegrals.

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

  unsigned int* indices = paramIndices + parameters;   
  devcomplex<fptype> ret = device_DalitzPlot_calcIntegrals(binCenterM12, binCenterM13, resonance_i, resonance_j, cudaArray, indices); 

  fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an event-weighted fit. 
  fakeEvt[indices[indices[0] + 2 + 0]] = binCenterM12;
  fakeEvt[indices[indices[0] + 2 + 1]] = binCenterM13;
  unsigned int numResonances = indices[2]; 
  int effFunctionIdx = parIndexFromResIndex_DP(numResonances); 
  fptype eff = callFunction(fakeEvt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 

  // Multiplication by eff, not sqrt(eff), is correct:
  // These complex numbers will not be squared when they
  // go into the integrals. They've been squared already,
  // as it were. 
  ret *= eff;
  return ret; 
}

SpecialResonanceCalculator::SpecialResonanceCalculator (int pIdx, unsigned int res_idx) 
  : resonance_i(res_idx)
  , parameters(pIdx)
{}

EXEC_TARGET devcomplex<fptype> SpecialResonanceCalculator::operator () (thrust::tuple<int, fptype*, int> t) const {
  // Calculates the BW values for a specific resonance. 
  devcomplex<fptype> ret;
  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t)); 

  unsigned int* indices = paramIndices + parameters;   // Jump to DALITZPLOT position within parameters array
  fptype m12 = evt[indices[2 + indices[0]]]; 
  fptype m13 = evt[indices[3 + indices[0]]];

  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3];  
  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 

  int parameter_i = parIndexFromResIndex_DP(resonance_i); // Find position of this resonance relative to DALITZPLOT start 

  unsigned int functn_i = indices[parameter_i+2];
  unsigned int params_i = indices[parameter_i+3];

  ret = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);
  //printf("Amplitude %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag); 
  return ret;
}

