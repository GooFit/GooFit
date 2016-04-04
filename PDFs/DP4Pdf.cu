#include "DP4Pdf.hh"
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

  int evtNum = (int) FLOOR(0.5 + evt[indices[7 + indices[0]]]); 
  // printf("%i\n",evtNum );
  devcomplex<fptype> totalAmp(0, 0);
  unsigned int cacheToUse    = indices[2]; 
  unsigned int numAmps       = indices[5]; 

  for (int i = 0; i < numAmps; ++i) {
    fptype amp_real = p[indices[6 + 2*i]];
    fptype amp_imag = p[indices[7 + 2*i]];

    devcomplex<fptype> matrixelement((Amps_DP[cacheToUse][evtNum*numAmps + i]).real,
             (Amps_DP[cacheToUse][evtNum*numAmps + i]).imag); 
    // printf("ddp %f, %f, %f, %f,\n",amp_real, amp_imag, matrixelement.real, matrixelement.imag );
    matrixelement.multiply(amp_real, amp_imag); 
    totalAmp += matrixelement;
  } 
  fptype ret = norm2(totalAmp); 
  int effFunctionIdx = 6 + 2*indices[3] + 2*indices[4] + 2*indices[5]; 
  fptype eff = callFunction(evt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  ret *= eff;

  // printf("test\n");
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_DP = device_DP; 

__host__ DPPdf::DPPdf (std::string n, 
                 std::vector<Variable*> observables,
                 // Variable* eventNumber, 
                 DecayInfo_DP* decay, 
                 GooPdf* efficiency)
  : GooPdf(0, n) 
  , decayInfo(decay)
  , _observables(observables)
  , cachedAMPs(0)
  , cachedResSF(0) 
  , forceRedoIntegrals(true)
  , totalEventSize(1 + observables.size()) // number of observables plus eventnumber
  , cacheToUse(0) 
  , SpinsCalculated(false)
  , devNormArray(0)
{
  for (std::vector<Variable*>::iterator obsIT = observables.begin(); obsIT != observables.end(); ++obsIT) {
    registerObservable(*obsIT);
  }
  // registerObservable(eventNumber); 

  std::vector<fptype>decayConstants;
  decayConstants.push_back(decayInfo->meson_radius);
  for(std::vector<fptype>::iterator pmIT = decayInfo->particle_masses.begin(); pmIT != decayInfo->particle_masses.end(); ++pmIT) {
    decayConstants.push_back(*pmIT);
  }
  
  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(decayConstants.size())); 
  MEMCPY_TO_SYMBOL(functorConstants, &decayConstants[0], decayConstants.size()*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);  
  static int cacheCount = 0; 
  cacheToUse = cacheCount++; 
  pindices.push_back(cacheToUse); 
  pindices.push_back(0); //#LS
  pindices.push_back(0); //#SF
  pindices.push_back(0); //#AMP 


  std::vector<unsigned int> ampidx;
  std::vector<unsigned int> nPermVec;
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
    nPermVec.push_back(decayInfo->amplitudes[i]->_nPerm);
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
  MEMCPY_TO_SYMBOL(AmpIndices, &(ampidxstart[0]), ampidxstart.size()*sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
  MEMCPY_TO_SYMBOL(AmpIndices, &(ampidx[0]), ampidx.size()*sizeof(unsigned int), ampidxstart.size()*sizeof(unsigned int), cudaMemcpyHostToDevice);
  
  pindices[2] = components.size(); 
  pindices[3] = SpinFactors.size();
  pindices[4] = AmpMap.size();

  for(int i = 0; i<components.size(); i++){
    reinterpret_cast<Lineshape*>(components[i])->setConstantIndex(cIndex);
    pindices.push_back(reinterpret_cast<Lineshape*>(components[i])->getFunctionIndex());
    pindices.push_back(reinterpret_cast<Lineshape*>(components[i])->getParameterIndex());

  }
  for (int i = 0; i < SpinFactors.size(); ++i)
  {
    pindices.push_back(SpinFactors[i]->getFunctionIndex());
    pindices.push_back(SpinFactors[i]->getParameterIndex());
    SpinFactors[i]->setConstantIndex(cIndex);
  }


  pindices.push_back(efficiency->getFunctionIndex());
  pindices.push_back(efficiency->getParameterIndex());
  components.push_back(efficiency); 


  GET_FUNCTION_ADDR(ptr_to_DP);
  initialise(pindices);

  redoIntegral = new bool[components.size() - 1];
  cachedMasses = new fptype[components.size() - 1];
  cachedWidths = new fptype[components.size() - 1];
  // lscalculators  = new LSCalculator*[components.size() - 1];
  // sfcalculators  = new SFCalculator*[SpinFactors.size()];

  for (int i = 0; i < components.size() - 1; ++i) {
    redoIntegral[i] = true;
    cachedMasses[i] = -1;
    cachedWidths[i] = -1; 
    lscalculators.push_back(new LSCalculator(parameters, i)); 
    
  }

  for (int i = 0; i < SpinFactors.size(); ++i) {
    sfcalculators.push_back(new SFCalculator(parameters, i)); 
    
  }
  for (int i = 0; i < AmpMap.size(); ++i)
  {
    integrators.push_back(new NormIntegrator(parameters, nPermVec[i]));
    AmpCalcs.push_back(new AmpCalc(ampidxstart[i], parameters, nPermVec[i]));
  }

  //printf("%i\n", parameters );

  addSpecialMask(PdfBase::ForceSeparateNorm); 
}

__host__ void DPPdf::setDataSize (unsigned int dataSize, unsigned int evtSize) {
  // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum = 6
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
  

  if (!devNormArray) {
    fptype* host_norm_array = new fptype[MCevents * (5 + 2*(components.size() - 1)  + SpinFactors.size() )];
    unsigned int offset = 5 + 2*(components.size() - 1)  + SpinFactors.size() ;
    // printf("devarray %i, %i\n",MCevents, offset );
    for (unsigned int i = 0; i < MCevents; ++i)
    {
      //m12 m34 cos12 cos34 phi weight
      host_norm_array[i*offset + 0] = hostphsp[i*5];    
      host_norm_array[i*offset + 1] = hostphsp[1+ i*5];
      host_norm_array[i*offset + 2] = hostphsp[2 + i*5];
      host_norm_array[i*offset + 3] = hostphsp[3 + i*5];
      host_norm_array[i*offset + 4] = hostphsp[4 + i*5];
      //host_norm_array[i*offset + 5] = hostphsp[5 + i*6];
    }
    // printf("%4g, %4g, %4g, %4g, %4g, %4g\n", host_norm_array[offset*(MCevents-1)], host_norm_array[offset*(MCevents-1)+1], host_norm_array[offset*(MCevents-1)+2], host_norm_array[offset*(MCevents-1)+3], host_norm_array[offset*(MCevents-1)+4], host_norm_array[offset*(MCevents-1)+5] );

    gooMalloc((void**) &devNormArray,( MCevents*offset*sizeof(fptype)) );
    MEMCPY(devNormArray, host_norm_array, (MCevents*offset*sizeof(fptype)), cudaMemcpyHostToDevice);
    delete[] host_norm_array; 
    delete[] hostphsp;
  }


  for (unsigned int i = 0; i < components.size() - 1; ++i) {
      redoIntegral[i] = forceRedoIntegrals; 
      if (!(components[i]->parametersChanged())) continue;
      redoIntegral[i] = true; 
      components[i]->storeParameters();
  }
  forceRedoIntegrals = false; 

  thrust::constant_iterator<fptype*> normaddress(devNormArray); 
  thrust::counting_iterator<int> binIndex(0); 
  thrust::constant_iterator<fptype*> dataArray(dev_event_array); 
  thrust::constant_iterator<int> eventSize(totalEventSize);
  thrust::counting_iterator<int> eventIndex(0); 

  if(!SpinsCalculated){
    for (int i = 0; i < SpinFactors.size() ; ++i) {
      unsigned int offset = components.size() -1;
      thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries , dataArray, eventSize)),
      strided_range<DEVICE_VECTOR<devcomplex<fptype> >::iterator>(cachedResSF->begin() + offset + i, 
                          cachedResSF->end(), 
                          (components.size() + SpinFactors.size() - 1)).begin(),
                          *(sfcalculators[i]));
    

    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(binIndex, normaddress)),
            thrust::make_zip_iterator(thrust::make_tuple(binIndex + MCevents, normaddress)),
            NormSpinCalculator(parameters, i));

    }
     SpinsCalculated = true;
  }

  for (int i = 0; i < components.size() -1 ; ++i) {
    if (redoIntegral[i]) {
      thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries , dataArray, eventSize)),
      strided_range<DEVICE_VECTOR<devcomplex<fptype> >::iterator>(cachedResSF->begin() + i, 
                        cachedResSF->end(), 
                        (components.size() + SpinFactors.size() - 1)).begin(), 
      *(lscalculators[i]));
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
      }
    }
    

  for (int i = 0; i < components.size() -1 ; ++i) {
      if(!redoIntegral[i]) continue;
      thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(binIndex, normaddress)),
              thrust::make_zip_iterator(thrust::make_tuple(binIndex + MCevents, normaddress)),
              NormLSCalculator(parameters, i));
    }

  fptype sumIntegral = 0;
  for(unsigned int i = 0; i<AmpCalcs.size(); ++i){
  sumIntegral += thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(binIndex, normaddress)),
          thrust::make_zip_iterator(thrust::make_tuple(binIndex + MCevents, normaddress)),
          *(integrators[i]),
          0.,
          thrust::plus<fptype>());
  }

  sumIntegral/=MCevents;
  host_normalisation[parameters] = 1.0/sumIntegral;
  // printf("end of normalise %f\n", sumIntegral);
  return sumIntegral;   
}

SFCalculator::SFCalculator (int pIdx, unsigned int sf_idx) 
  : _spinfactor_i(sf_idx)
  , _parameters(pIdx)
{}

EXEC_TARGET devcomplex<fptype> SFCalculator::operator () (thrust::tuple<int, fptype*, int> t) const {
  
  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t)); 

  unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
  int parameter_i = 6 + (2 * indices[5]) + (indices[3] * 2) + (_spinfactor_i * 2) ; // Find position of this resonance relative to DALITZPLOT start 
  unsigned int functn_i = indices[parameter_i];
  unsigned int params_i = indices[parameter_i+1];
  
  fptype m12 = evt[indices[2 + indices[0]]]; 
  fptype m34 = evt[indices[3 + indices[0]]];
  fptype cos12 = evt[indices[4 + indices[0]]];
  fptype cos34 = evt[indices[5 + indices[0]]];
  fptype phi = evt[indices[6 + indices[0]]];

  fptype vecs[16];
  get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi); 
  // printf("%f, %f, %f, %f, %f \n",m12, m34, cos12, cos34, phi );
  // printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
  // printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
  // printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
  // printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);

  spin_function_ptr func = reinterpret_cast<spin_function_ptr>(device_function_table[functn_i]);
  fptype sf = (*func)(vecs, paramIndices+params_i);
  // printf("SpinFactors: %f\n", sf );
  return devcomplex<fptype>(sf, 0);
}

NormSpinCalculator::NormSpinCalculator (int pIdx, unsigned int sf_idx) 
  : _spinfactor_i(sf_idx)
  , _parameters(pIdx)
{}

EXEC_TARGET fptype NormSpinCalculator::operator () (thrust::tuple<int, fptype*> t) const {

  unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
  unsigned int numLS    = indices[3];
  unsigned int numSF    = indices[4];
  unsigned int numAmps  = indices[5];
  int parameter_i = 6 + (2 * numAmps) + (numLS * 2) + (_spinfactor_i * 2) ; // Find position of this resonance relative to DALITZPLOT start 
  unsigned int offset = 5 + 2 * numLS + numSF;
  unsigned int functn_i = indices[parameter_i];
  unsigned int params_i = indices[parameter_i+1];
  
  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * offset); 
  fptype m12    = evt[0]; 
  fptype m34    = evt[1];
  fptype cos12  = evt[2];
  fptype cos34  = evt[3];
  fptype phi    = evt[4];

  fptype vecs[16];
  get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi); 

//   printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
//   printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
//   printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
//   printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);
// }
  spin_function_ptr func = reinterpret_cast<spin_function_ptr>(device_function_table[functn_i]);
  fptype sf = (*func)(vecs, paramIndices+params_i);
  fptype* spinfactor = evt + 5 + 2*numLS + _spinfactor_i;
  spinfactor[0] = sf;
  // printf("%f\n",sf );

  THREAD_SYNCH
  return sf;
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
  int parameter_i = 6 + (2 * indices[5]) + (_resonance_i * 2) ; // Find position of this resonance relative to DALITZPLOT start 
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

  if (pair < 2){
    fptype mres = pair==0 ? m12 : m34;
    fptype d1 = pair==0 ? m1 : m3;
    fptype d2 = pair==0 ? m2 : m4;
    ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
    // printf("LS %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
  }
  else{ 
    fptype vecs[16];
    get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi); 
    fptype d1,d2;  
    fptype mres = getmass(pair, d1 , d2, vecs, m1, m2, m3, m4);
    ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
  }

  //if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  //printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag); 
  //printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
  return ret;
}

NormLSCalculator::NormLSCalculator (int pIdx, unsigned int res_idx) 
  : _resonance_i(res_idx)
  , _parameters(pIdx)
{}

EXEC_TARGET int NormLSCalculator::operator () (thrust::tuple<int, fptype*> t) const {
  // Calculates the BW values for a specific resonance. 
  devcomplex<fptype> ret;
  
  unsigned int* indices = paramIndices + _parameters;   // Jump to DALITZPLOT position within parameters array
  unsigned int numLS    = indices[3];
  unsigned int numSF    = indices[4];
  unsigned int numAmps  = indices[5];
  int parameter_i = 6 + (2 * numAmps) + (_resonance_i * 2); // Find position of this resonance relative to DALITZPLOT start 
  unsigned int offset = 5 + 2 * numLS + numSF;
  unsigned int functn_i = indices[parameter_i];
  unsigned int params_i = indices[parameter_i+1];
  unsigned int pair = (paramIndices+params_i)[5];
  
  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * offset); 

  fptype m1  = functorConstants[indices[1] + 2]; 
  fptype m2  = functorConstants[indices[1] + 3]; 
  fptype m3  = functorConstants[indices[1] + 4];  
  fptype m4  = functorConstants[indices[1] + 5];
  
  fptype m12 = evt[0]; 
  fptype m34 = evt[1];
  fptype cos12 = evt[2];
  fptype cos34 = evt[3];
  fptype phi = evt[4];

  if (pair < 2){
    fptype mres = pair==0 ? m12 : m34;
    fptype d1 = pair==0 ? m1 : m3;
    fptype d2 = pair==0 ? m2 : m4;
    ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
  }
  else{ 
    fptype vecs[16];
    get4Vecs(vecs, indices[1], m12, m34, cos12, cos34, phi); 
    fptype d1,d2;  
    fptype mres = getmass(pair, d1 , d2, vecs, m1, m2, m3, m4);
    ret = getResonanceAmplitude(mres, d1, d2, functn_i, params_i);
  }
  // printf("%f, %f, %f, %f, %f \n",m12, m34, cos12, cos34, phi );
  // printf("%i, %i, %i, %i, %i \n",numLS, numSF, numAmps, offset, evtNum );
  // printf("NLS %f, %f\n",ret.real, ret.imag);
  evt[5 + _resonance_i * 2] = ret.real;
  evt[6 + _resonance_i * 2] = ret.imag;

  //if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  //printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag); 
  //printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
  THREAD_SYNCH
  return 0;
}

AmpCalc::AmpCalc(unsigned int AmpIdx, unsigned int pIdx, unsigned int nPerm)
  : _nPerm(nPerm)
  , _AmpIdx(AmpIdx)
  , _parameters(pIdx)
  {}


EXEC_TARGET devcomplex<fptype> AmpCalc::operator() (thrust::tuple<int, fptype*, int> t) const {
  unsigned int * indices = paramIndices + _parameters;
  unsigned int cacheToUse = indices[2];
  unsigned int totalLS = indices[3];
  unsigned int totalSF = indices[4];
  unsigned int totalAMP = indices[5];
  unsigned int offset = totalLS + totalSF;
  unsigned int numLS = AmpIndices[totalAMP + _AmpIdx];
  unsigned int numSF = AmpIndices[totalAMP + _AmpIdx + 1];
  unsigned int evtNum = thrust::get<0>(t);

  devcomplex<fptype> returnVal(0,0);
  unsigned int SF_step = numSF/_nPerm;
  unsigned int LS_step = numLS/_nPerm;

  for (int i = 0; i < _nPerm; ++i)
  {
    devcomplex<fptype> ret(1,0);
    for (int j = i*LS_step; j < (i+1)*LS_step; ++j){
      ret *= (cResSF[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 2 + j]]);
    }
    for (int j = i*SF_step; j < (i+1)*SF_step; ++j){
      ret *= (cResSF[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 2 + numLS + j]].real);
    }
    returnVal += ret;
  }
 
  return (1/SQRT((fptype)(_nPerm))) * returnVal;
}

NormIntegrator::NormIntegrator(unsigned int pIdx, unsigned int nPerm)
  : _parameters(pIdx)
  , _nPerm(nPerm)
  {}


EXEC_TARGET fptype NormIntegrator::operator() (thrust::tuple<int, fptype*> t) const {
  unsigned int * indices = paramIndices + _parameters;
  unsigned int totalLS = indices[3];
  unsigned int totalSF = indices[4];
  unsigned int totalAMP = indices[5];
  unsigned int offset = 5 + 2 * totalLS + totalSF;

  unsigned int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * offset); 

  devcomplex<fptype> returnVal(0,0);

  for (int amp = 0; amp < totalAMP; ++amp)
  {
    unsigned int ampidx =  AmpIndices[amp];
    unsigned int numLS = AmpIndices[totalAMP + ampidx];
    unsigned int numSF = AmpIndices[totalAMP + ampidx + 1];
    unsigned int SF_step = numSF/_nPerm;
    unsigned int LS_step = numLS/_nPerm;
    devcomplex<fptype> ret2(0,0);

    for (int j = 0; j < _nPerm; ++j){  
      devcomplex<fptype> ret(1,0);
      for (int i = j*LS_step; i < (j+1)*LS_step; ++i){
        devcomplex<fptype> matrixelement(evt[5 + 2 * AmpIndices[totalAMP + ampidx + 2 + i]],evt[5 + 2 * AmpIndices[totalAMP + ampidx + 2 + i] + 1]); 
        ret *= matrixelement; 
      }
      for (int i = j*SF_step; i < (j+1)*SF_step; ++i){
        devcomplex<fptype> matrixelement(evt[5 + 2 * totalLS + AmpIndices[totalAMP + ampidx + 2 + numLS + i]],0); 
        ret *= matrixelement; 
      }
      ret2 += ret;
    }

    ret2 *= (1/SQRT((fptype)(_nPerm)));
    fptype amp_real = cudaArray[indices[2*amp + 5]];
    fptype amp_imag = cudaArray[indices[2*amp + 6]];
    ret2.multiply(amp_real, amp_imag); 
    returnVal += ret2;
  }
  return norm2(returnVal);
}




