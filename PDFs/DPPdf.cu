#include "DPPdf.hh"
#include <complex>
using std::complex; 

//const int resonanceOffset_DP = 4; // Offset of the first resonance into the parameter index array 
// Offset is number of parameters, constant index, number of resonances (not calculable 
// from nP because we don't know what the efficiency might need), and cache index. Efficiency 
// parameters are after the resonance information. 

// The function of this array is to hold all the cached waves; specific 
// waves are recalculated when the corresponding resonance mass or width 
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone! 
MEM_DEVICE devcomplex<fptype>* cResSF[10]; 
MEM_DEVICE devcomplex<fptype>* Amps_DP[10]; 
MEM_CONSTANT unsigned int AmpIndices[100];

EXEC_TARGET fptype device_DP (fptype* evt, fptype* p, unsigned int* indices) {
  //printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real, totalAmp.imag); 

  // printf("test\n");
  return 1.0; 
}

MEM_DEVICE device_function_ptr ptr_to_DP = device_DP; 

__host__ DPPdf::DPPdf (std::string n, 
							   std::vector<Variable*> observables,
							   Variable* eventNumber, 
							   DecayInfo_DP* decay, 
							   GooPdf* efficiency)
  : GooPdf(0, n) 
  , decayInfo(decay)
  , _observables(observables)
  , dalitzNormRange(0)
  , cachedAMPs(0)
  , cachedResSF(0) 
  , integrals(0)
  , forceRedoIntegrals(true)
  , totalEventSize(1 + observables.size()) // number of observables plus eventnumber
  , cacheToUse(0) 
  , integrators(0)
  , calculators(0) 
  , SpinsCalculated(false)
{
  for (std::vector<Variable*>::iterator obsIT = observables.begin(); obsIT != observables.end(); ++obsIT) {
    registerObservable(*obsIT);
  }
  registerObservable(eventNumber); 

  fptype decayConstants[1 + decayInfo->particle_masses.size()];
  int a = 0;
  decayConstants[a] = decayInfo->meson_radius;
  for(std::map<unsigned int, fptype>::iterator pmIT = decayInfo->particle_masses.begin(); pmIT != decayInfo->particle_masses.end(); ++pmIT) {
    massmap[pmIT->first] = a + 1;
    decayConstants[a + 1] = pmIT->second;
    a++;
  }
  
  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(a)); 
  MEMCPY_TO_SYMBOL(functorConstants, decayConstants, a*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);  
  static int cacheCount = 0; 
  cacheToUse = cacheCount++; 
  pindices.push_back(cacheToUse); 
  pindices.push_back(0); //#LS
  pindices.push_back(0); //#SF
  pindices.push_back(0); //#AMP 


  std::vector<unsigned int> ampidx;
  std::vector<unsigned int> ampidxstart;
  for(int i=0; i<decayInfo->amplitudes.size(); i++){
    AmpMap[decayInfo->amplitudes[i]->_uniqueDecayStr] =  std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));
    for(std::map<std::string, Lineshape*>::iterator LSIT = decayInfo->amplitudes[i]->_LS.begin(); LSIT != decayInfo->amplitudes[i]->_LS.end(); ++LSIT) {
      if(compMap.find(LSIT->first) != compMap.end()){
        AmpMap[decayInfo->amplitudes[i]->_uniqueDecayStr].first.push_back(compMap.find(LSIT->first)->second);
      }else{
        components.push_back(LSIT->second);
        compMap[LSIT->first] = components.size() - 1; 
        AmpMap[decayInfo->amplitudes[i]->_uniqueDecayStr].first.push_back(compMap.find(LSIT->first)->second);
      }
    }
    for(std::map<std::string, SpinFactor*>::iterator SFIT = decayInfo->amplitudes[i]->_SF.begin(); SFIT != decayInfo->amplitudes[i]->_SF.end(); ++SFIT) {
      if(SpinMap.find(SFIT->first) != SpinMap.end()){
        AmpMap[decayInfo->amplitudes[i]->_uniqueDecayStr].second.push_back(SpinMap.find(SFIT->first)->second);
      }else{
        SpinFactors.push_back(SFIT->second);
        SpinMap[SFIT->first] = SpinFactors.size() - 1; 
        AmpMap[decayInfo->amplitudes[i]->_uniqueDecayStr].second.push_back(SpinMap.find(SFIT->first)->second);
      }
    }

    pindices.push_back(registerParameter(decayInfo->amplitudes[i]->_ar));
    pindices.push_back(registerParameter(decayInfo->amplitudes[i]->_ai));
    
    // AmpCalcs.push_back(new AmpCalc(ampidx.size(), i));
    ampidxstart.push_back(ampidx.size());
    std::vector<unsigned int> ls = AmpMap[decayInfo->amplitudes[i]->_uniqueDecayStr].first;
    std::vector<unsigned int> sf = AmpMap[decayInfo->amplitudes[i]->_uniqueDecayStr].second;
    ampidx.push_back(ls.size());
    ampidx.push_back(sf.size());
    ampidx.insert(ampidx.end(), ls.begin(), ls.end());
    ampidx.insert(ampidx.end(), sf.begin(), sf.end());
  }

  MEMCPY_TO_SYMBOL(AmpIndices, &(ampidx[0]), ampidx.size()*sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
  
  pindices[2] = components.size(); 
  pindices[3] = SpinFactors.size();
  pindices[4] = AmpMap.size();

  for(int i = 0; i<components.size(); i++){
    reinterpret_cast<Lineshape*>(components[i])->setMotherIndex(massmap[reinterpret_cast<Lineshape*>(components[i])->_mother_pdg]);
    reinterpret_cast<Lineshape*>(components[i])->setConstantIndex(cIndex);
    pindices.push_back(reinterpret_cast<Lineshape*>(components[i])->getFunctionIndex());
    pindices.push_back(reinterpret_cast<Lineshape*>(components[i])->getParameterIndex());

  }
  for (int i = 0; i < SpinFactors.size(); ++i)
  {
    pindices.push_back(SpinFactors[i]->getFunctionIndex());
    pindices.push_back(SpinFactors[i]->getParameterIndex());
  }


  pindices.push_back(efficiency->getFunctionIndex());
  pindices.push_back(efficiency->getParameterIndex());
  components.push_back(efficiency); 


  GET_FUNCTION_ADDR(ptr_to_DP);
  initialise(pindices);

  redoIntegral = new bool[components.size() - 1];
  cachedMasses = new fptype[components.size() - 1];
  cachedWidths = new fptype[components.size() - 1];
  integrals    = new devcomplex<fptype>**[AmpMap.size()];
  integrators  = new SpecialIntegrator**[AmpMap.size()];
  calculators  = new LSCalculator*[components.size() - 1];

  for (int i = 0; i < components.size() - 1; ++i) {
    redoIntegral[i] = true;
    cachedMasses[i] = -1;
    cachedWidths[i] = -1; 
    calculators[i]  = new LSCalculator(parameters, i); 
    
  }

  for (int i = 0; i < AmpMap.size(); ++i)
  {
    integrators[i]  = new SpecialIntegrator*[AmpMap.size()];
    integrals[i]    = new devcomplex<fptype>*[AmpMap.size()];
    AmpCalcs.push_back(new AmpCalc(ampidxstart[i], parameters, i));
    for (int j = 0; j < AmpMap.size(); ++j) {
      integrals[i][j]   = new devcomplex<fptype>(0, 0); 
      integrators[i][j] = new SpecialIntegrator(parameters, i, j, ampidxstart[i], ampidxstart[j]); 
    }
  }

  // for (int i = 0; i < SpinFactors.size(); ++i)
  // {
  //   SpinFactors[i]->SetParameterIdx(parameters);
  //   SpinFactors[i]->resolveMassIdx(massmap);
  // }
  printf("%i\n", parameters );

  addSpecialMask(PdfBase::ForceSeparateNorm); 
}

__host__ void DPPdf::setDataSize (unsigned int dataSize, unsigned int evtSize) {
  // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum
  totalEventSize = evtSize;
  assert(totalEventSize >= 3); 

  if (cachedResSF) delete cachedResSF;
  if (cachedAMPs) delete cachedAMPs;

  numEntries = dataSize; 
  cachedResSF = new DEVICE_VECTOR<devcomplex<fptype> >(dataSize*(components.size() + SpinFactors.size() - 1)); //   -1 because 1 component is efficiency
  void* dummy = thrust::raw_pointer_cast(cachedResSF->data()); 
  MEMCPY_TO_SYMBOL(cResSF, &dummy, sizeof(devcomplex<fptype>*), cacheToUse*sizeof(devcomplex<fptype>*), cudaMemcpyHostToDevice); 
  
  cachedAMPs = new DEVICE_VECTOR<devcomplex<fptype> >(dataSize*(AmpCalcs.size())); 
  void* dummy2 = thrust::raw_pointer_cast(cachedAMPs->data()); 
  MEMCPY_TO_SYMBOL(Amps_DP, &dummy2, sizeof(devcomplex<fptype>*), cacheToUse*sizeof(devcomplex<fptype>*), cudaMemcpyHostToDevice); 

  setForceIntegrals(); 
}

__host__ fptype DPPdf::normalise () const {
  recursiveSetNormalisation(1); // Not going to normalise efficiency, 
  // so set normalisation factor to 1 so it doesn't get multiplied by zero. 
  // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency, 
  // don't get zeroes through multiplying by the normFactor. 
  MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  
  int totalBins = 1;
  for (int i = 0; i < _observables.size(); ++i){totalBins *= _observables[i]->numbins;}

  if (!dalitzNormRange) {
    gooMalloc((void**) &dalitzNormRange, (_observables.size() * 3)*sizeof(fptype));
    fptype* host_norms = new fptype[3 * _observables.size()];
    for (int i = 0; i < _observables.size(); ++i){
      host_norms[3*i] = _observables[i]->lowerlimit;
      host_norms[3*i+1] = _observables[i]->upperlimit;
      host_norms[3*i+2] = _observables[i]->numbins;
    }
    MEMCPY(dalitzNormRange, host_norms, (_observables.size() * 3)*sizeof(fptype), cudaMemcpyHostToDevice);
    delete[] host_norms; 
  }

  for (unsigned int i = 0; i < components.size() - 1; ++i) {
      redoIntegral[i] = forceRedoIntegrals; 
      if (!(components[i]->parametersChanged())) continue;
      redoIntegral[i] = true; 
      components[i]->storeParameters();
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

  // if(!SpinsCalculated){
  //   for (int i = 0; i < SpinFactors.size() ; ++i) {
  //     unsigned int offset = components.size() -1;
  //     thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
  //     thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries , dataArray, eventSize)),
  //     strided_range<DEVICE_VECTOR<devcomplex<fptype> >::iterator>(cachedResSF->begin() + offset + i, 
  //                         cachedResSF->end(), 
  //                         (components.size() + SpinFactors.size() - 1)).begin(),
  //                         *(SpinFactors[i]));
  //   }
  //    SpinsCalculated = true;
  // }

  for (int i = 0; i < components.size() -1 ; ++i) {
    if (redoIntegral[i]) {
      thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries , dataArray, eventSize)),
      strided_range<DEVICE_VECTOR<devcomplex<fptype> >::iterator>(cachedResSF->begin() + i, 
                        cachedResSF->end(), 
                        (components.size() + SpinFactors.size() - 1)).begin(), 
      *(calculators[i]));
    }
  }

  std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int> > >::const_iterator AmpMapIt = AmpMap.begin();
  for (int i = 0; i < AmpCalcs.size(); ++i) {
    std::vector<unsigned int> redoidx((*AmpMapIt).second.first);
    bool redo = false;
    for (int j = 0; j < redoidx.size(); ++j){
      if(!redoIntegral[redoidx[j]]) continue;
      redo = true;
      break;
    }
    if (redo) {
      thrust::transform(eventIndex, eventIndex + numEntries,
                        strided_range<DEVICE_VECTOR<devcomplex<fptype> >::iterator>(cachedAMPs->begin() + i, 
                        cachedAMPs->end(), AmpCalcs.size()).begin(), 
                        *(AmpCalcs[i]));


      for (int j = 0; j < AmpCalcs.size(); ++j)
      {
        devcomplex<fptype> dummy(0, 0);
        thrust::plus<devcomplex<fptype> > complexSum; 
        (*(integrals[i][j])) = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
                    thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
                    *(integrators[i][j]), 
                    dummy, 
                    complexSum); 
        // printf("inside of integration loop %i %f \n", totalBins, integrals[i][j]->imag );
        }
      }
    }
  complex<fptype> sumIntegral(0, 0);
  for (unsigned int i = 0; i < AmpCalcs.size(); ++i) {
    int param_i = parameters + 6 + 2 * i; 
    complex<fptype> amplitude_i(host_params[host_indices[param_i]], host_params[host_indices[param_i + 1]]);
    for (unsigned int j = 0; j < AmpCalcs.size(); ++j) {
      int param_j = parameters + 6 + 2 * j; 
      complex<fptype> amplitude_j(host_params[host_indices[param_j]], -host_params[host_indices[param_j + 1]]); 
      // Notice complex conjugation

      printf("%f %f %f %f %f %f\n", amplitude_i.real(), amplitude_i.imag(), amplitude_j.real(), amplitude_j.imag(), (*(integrals[i][j])).real, (*(integrals[i][j])).imag );
      sumIntegral += (amplitude_i * amplitude_j * complex<fptype>((*(integrals[i][j])).real, (*(integrals[i][j])).imag)); 
    }
  }

  fptype ret = real(sumIntegral); // That complex number is a square, so it's fully real
  double binSizeFactor = 1;
  for (int i = 0; i < _observables.size(); ++i){
    binSizeFactor *= (_observables[i]->upperlimit-_observables[i]->lowerlimit) / _observables[i]->numbins;
  }
  ret *= binSizeFactor;

  host_normalisation[parameters] = 1.0/ret;
  printf("%f, %f \n", ret, binSizeFactor);
  return ret;   
}

SFCalculator::SFCalculator (int pIdx, unsigned int sf_idx) 
  : _spinfactor_i(sf_idx)
  , _parameters(pIdx)
{}

EXEC_TARGET devcomplex<fptype> SFCalculator::operator () (thrust::tuple<int, fptype*, int> t) const {
  // Calculates the BW values for a specific resonance. 
  devcomplex<fptype> ret(1,0);
  
  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t)); 

  unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
  fptype m12 = evt[indices[2 + indices[0]]]; 
  fptype m13 = evt[indices[3 + indices[0]]];
  fptype motherMass = functorConstants[indices[1] + 1]; 
  fptype daug1Mass  = functorConstants[indices[1] + 2]; 
  fptype daug2Mass  = functorConstants[indices[1] + 3]; 
  fptype daug3Mass  = functorConstants[indices[1] + 4];  
  //if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 

  int parameter_i = 6 + (2 * indices[5]) + (_resonance_i * 2) ; // Find position of this resonance relative to DALITZPLOT start 

  unsigned int functn_i = indices[parameter_i];
  unsigned int params_i = indices[parameter_i+1];

  //printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag); 
  //printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
  return ret;
}


LSCalculator::LSCalculator (int pIdx, unsigned int res_idx) 
  : _resonance_i(res_idx)
  , _parameters(pIdx)
{}

EXEC_TARGET devcomplex<fptype> LSCalculator::operator () (thrust::tuple<int, fptype*, int> t) const {
  // Calculates the BW values for a specific resonance. 
  devcomplex<fptype> ret;
  
  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t)); 

  unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
  fptype m12 = evt[indices[2 + indices[0]]]; 
  fptype m13 = evt[indices[3 + indices[0]]];
  fptype motherMass = functorConstants[indices[1] + 1]; 
  fptype daug1Mass  = functorConstants[indices[1] + 2]; 
  fptype daug2Mass  = functorConstants[indices[1] + 3]; 
  fptype daug3Mass  = functorConstants[indices[1] + 4];  
  //if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 

  int parameter_i = 6 + (2 * indices[5]) + (_resonance_i * 2) ; // Find position of this resonance relative to DALITZPLOT start 

  unsigned int functn_i = indices[parameter_i];
  unsigned int params_i = indices[parameter_i+1];

  ret = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);

  //printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag); 
  //printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
  return ret;
}

AmpCalc::AmpCalc(unsigned int AmpIdx, unsigned int pIdx, unsigned int coeffIdx)
  : _coeffIdx(coeffIdx)
  , _AmpIdx(AmpIdx)
  , _parameters(pIdx)
  {}


EXEC_TARGET devcomplex<fptype> AmpCalc::operator() (thrust::tuple<int, fptype*, int> t) const {
  devcomplex<fptype> ret(0,0);
  unsigned int * indices = paramIndices + _parameters;
  unsigned int cacheToUse = indices[2];
  unsigned int totalLS = indices[3];
  unsigned int totalSF = indices[4];
  unsigned int offset = totalLS + totalSF;
  unsigned int numLS = AmpIndices[_AmpIdx];
  unsigned int numSF =AmpIndices[_AmpIdx + 1];
  unsigned int evtNum = thrust::get<0>(t);

  for (int i = 0; i < (numLS + numSF); ++i)
  {
    unsigned int a = (i<numLS ? 0 : totalLS);
    fptype amp_real = cudaArray[indices[2*_coeffIdx + 6]];
    fptype amp_imag = cudaArray[indices[2*_coeffIdx + 7]];

    devcomplex<fptype> matrixelement((cResSF[cacheToUse][evtNum*offset + a + AmpIndices[_AmpIdx + 2 + i]]).real,
             (cResSF[cacheToUse][evtNum*offset + a + AmpIndices[_AmpIdx + 2 + i]]).imag); 
    matrixelement.multiply(amp_real, amp_imag); 
    ret *= matrixelement; 
  }
  //printf("test %i, %i, %i, %i, %i, %i, %i, \n", cacheToUse, totalLS, totalSF, offset, numLS, numSF, evtNum);
  return ret;
}

SpecialIntegrator::SpecialIntegrator (int pIdx, unsigned int ri, unsigned int rj, unsigned int starti, unsigned int startj)
  : _parameters(pIdx)
  , _amp_i(ri)
  , _amp_j(rj)
  , _starti(starti)
  , _startj(startj)
{}

EXEC_TARGET devcomplex<fptype> SpecialIntegrator::operator () (thrust::tuple<int, fptype*> t) const{
  //unsigned int NumObs = indices[indices[0] + 1];
  int globalBinNumber  = thrust::get<0>(t);
  fptype lowerBoundM12 = thrust::get<1>(t)[0];
  fptype upperBoundM12 = thrust::get<1>(t)[1]; 
  int numBinsM12       = (int) FLOOR(thrust::get<1>(t)[2] + 0.5); 
  printf("%i, %i, %f, %f\n", globalBinNumber, numBinsM12, lowerBoundM12, upperBoundM12);
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

  unsigned int* indices = paramIndices + _parameters;   
  fptype motherMass = functorConstants[indices[1] + 1]; 
  fptype daug1Mass  = functorConstants[indices[1] + 2]; 
  fptype daug2Mass  = functorConstants[indices[1] + 3]; 
  fptype daug3Mass  = functorConstants[indices[1] + 4];  

  // fptype binCenterM12 = 1.5;
  // fptype binCenterM13 = 1.6;
  

  if (!inDalitz(binCenterM12, binCenterM13, motherMass, daug1Mass, daug2Mass, daug3Mass)){
    printf("not accepted binm12 %f, binCenterM13 %f\n", binCenterM12, binCenterM13);
    return devcomplex<fptype>(0,0);
  }
  printf("! accepted binm12 %f, binCenterM13 %f\n", binCenterM12, binCenterM13);
  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - binCenterM12 - binCenterM13; 

  unsigned int numLSi = AmpIndices[_starti];
  unsigned int numSFi =AmpIndices[_starti + 1];
  unsigned int numLSj = AmpIndices[_startj];
  unsigned int numSFj =AmpIndices[_startj + 1];
  unsigned int offset = 6 + indices[5] * 2;

  // printf("%i, %i\n", numLSi, AmpIndices[0] );
  devcomplex<fptype> ret (1,0);
  for (int i = 0; i < numLSi; ++i)
  {
    unsigned int fcnt = offset + 2 * AmpIndices[_starti + 2 + i];
    ret *= getResonanceAmplitude(binCenterM12, binCenterM13, m23, indices[fcnt], indices[fcnt +1 ]);
  }
  for (int i = 0; i < numLSj; ++i)
  {
    unsigned int fcnt = offset + 2 * AmpIndices[_startj + 2 + i];
    ret *= conj(getResonanceAmplitude(binCenterM12, binCenterM13, m23, indices[fcnt], indices[fcnt +1 ]));
  }
  fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an event-weighted fit. 
  fakeEvt[indices[indices[0] + 2 + 0]] = binCenterM12;
  fakeEvt[indices[indices[0] + 2 + 1]] = binCenterM13;
  int effFunctionIdx = offset + (2*indices[3]) + (indices[4] * 2) ; 
  fptype eff = callFunction(fakeEvt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  ret *= eff;
  printf("ret %f %f %f %f %f\n",binCenterM12, binCenterM13, ret.real, ret.imag, eff );
  return ret;
}


