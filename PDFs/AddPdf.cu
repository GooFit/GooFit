#include "AddPdf.hh"

EXEC_TARGET fptype device_AddPdfs (fptype* evt, fptype* p, unsigned int* indices) { 
  int numParameters = indices[0]; 
  fptype ret = 0;
  fptype totalWeight = 0; 
  for (int i = 1; i < numParameters-3; i += 3) {
    totalWeight += p[indices[i+2]];
    fptype curr = callFunction(evt, indices[i], indices[i+1]); 
    fptype weight = p[indices[i+2]];
    ret += weight * curr * normalisationFactors[indices[i+1]]; 

    //if ((gpuDebug & 1) && (0 == threadIdx.x) && (0 == blockIdx.x)) 
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("Add comp %i: %f * %f * %f = %f (%f)\n", i, weight, curr, normalisationFactors[indices[i+1]], weight*curr*normalisationFactors[indices[i+1]], ret); 

  }
  // numParameters does not count itself. So the array structure for two functions is
  // nP | F P w | F P
  // in which nP = 5. Therefore the parameter index for the last function pointer is nP, and the function index is nP-1. 
  //fptype last = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[numParameters-1]])))(evt, p, paramIndices + indices[numParameters]);
  fptype last = callFunction(evt, indices[numParameters - 1], indices[numParameters]);
  ret += (1 - totalWeight) * last * normalisationFactors[indices[numParameters]]; 

  //if ((threadIdx.x < 50) && (isnan(ret))) printf("NaN final component %f %f\n", last, totalWeight); 

  //if ((gpuDebug & 1) && (0 == threadIdx.x) && (0 == blockIdx.x)) 
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //printf("Add final: %f * %f * %f = %f (%f)\n", (1 - totalWeight), last, normalisationFactors[indices[numParameters]], (1 - totalWeight) *last* normalisationFactors[indices[numParameters]], ret); 
  
  return ret; 
}

EXEC_TARGET fptype device_AddPdfsExt (fptype* evt, fptype* p, unsigned int* indices) { 
  // numParameters does not count itself. So the array structure for two functions is
  // nP | F P w | F P w
  // in which nP = 6. 

  int numParameters = indices[0]; 
  fptype ret = 0;
  fptype totalWeight = 0; 
  for (int i = 1; i < numParameters; i += 3) {    
    //fptype curr = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[i]])))(evt, p, paramIndices + indices[i+1]);
    fptype curr = callFunction(evt, indices[i], indices[i+1]); 
    fptype weight = p[indices[i+2]];
    ret += weight * curr * normalisationFactors[indices[i+1]]; 

    totalWeight += weight; 
    //if ((gpuDebug & 1) && (threadIdx.x == 0) && (0 == blockIdx.x)) 
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("AddExt: %i %E %f %f %f %f %f %f\n", i, curr, weight, ret, totalWeight, normalisationFactors[indices[i+1]], evt[0], evt[8]);
  }
  ret /= totalWeight; 
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //if ((gpuDebug & 1) && (threadIdx.x == 0) && (0 == blockIdx.x)) 
  //printf("AddExt result: %f\n", ret); 
  
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_AddPdfs = device_AddPdfs; 
MEM_DEVICE device_function_ptr ptr_to_AddPdfsExt = device_AddPdfsExt; 

AddPdf::AddPdf (std::string n, std::vector<Variable*> weights, std::vector<PdfBase*> comps) 
  : GooPdf(0, n) 
  , extended(true)
{

  assert((weights.size() == comps.size()) || (weights.size() + 1 == comps.size())); 

  // Indices stores (function index)(function parameter index)(weight index) triplet for each component. 
  // Last component has no weight index unless function is extended. 
  for (std::vector<PdfBase*>::iterator p = comps.begin(); p != comps.end(); ++p) {
    components.push_back(*p); 
    assert(components.back()); 
  }

  getObservables(observables); 

  std::vector<unsigned int> pindices;
  for (unsigned int w = 0; w < weights.size(); ++w) {
    assert(components[w]);
    pindices.push_back(components[w]->getFunctionIndex());
    pindices.push_back(components[w]->getParameterIndex());
    pindices.push_back(registerParameter(weights[w])); 
  }
  assert(components.back()); 
  if (weights.size() < components.size()) {
    pindices.push_back(components.back()->getFunctionIndex());
    pindices.push_back(components.back()->getParameterIndex());
    extended = false; 
  }


  if (extended) MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_AddPdfsExt, sizeof(void*), 0, cudaMemcpyDeviceToHost);
  else MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_AddPdfs, sizeof(void*), 0, cudaMemcpyDeviceToHost);

  initialise(pindices); 
} 


AddPdf::AddPdf (std::string n, Variable* frac1, PdfBase* func1, PdfBase* func2) 
  : GooPdf(0, n) 
  , extended(false)
{
  // Special-case constructor for common case of adding two functions.
  components.push_back(func1);
  components.push_back(func2);
  getObservables(observables); 

  std::vector<unsigned int> pindices;
  pindices.push_back(func1->getFunctionIndex());
  pindices.push_back(func1->getParameterIndex());
  pindices.push_back(registerParameter(frac1)); 

  pindices.push_back(func2->getFunctionIndex());
  pindices.push_back(func2->getParameterIndex());
    
  MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_AddPdfs, sizeof(void*), 0, cudaMemcpyDeviceToHost);

  initialise(pindices); 
} 

__host__ fptype AddPdf::normalise () const {
  //if (cpuDebug & 1) std::cout << "Normalising AddPdf " << getName() << std::endl;

  fptype ret = 0;
  fptype totalWeight = 0; 
  for (unsigned int i = 0; i < components.size()-1; ++i) {
    fptype weight = host_params[host_indices[parameters + 3*(i+1)]]; 
    totalWeight += weight;
    fptype curr = components[i]->normalise(); 
    ret += curr*weight;
  }
  fptype last = components.back()->normalise(); 
  if (extended) {
    fptype lastWeight = host_params[host_indices[parameters + 3*components.size()]];
    totalWeight += lastWeight;
    ret += last * lastWeight; 
    ret /= totalWeight; 
  }
  else {
    ret += (1 - totalWeight) * last;
  }
  host_normalisation[parameters] = 1.0; 

  if (getSpecialMask() & PdfBase::ForceCommonNorm) {
    // Want to normalise this as 
    // (f1 A + (1-f1) B) / int (f1 A + (1-f1) B) 
    // instead of default 
    // (f1 A / int A) + ((1-f1) B / int B).

    for (unsigned int i = 0; i < components.size(); ++i) {
      host_normalisation[components[i]->getParameterIndex()] = (1.0 / ret);
    }
  }

  //if (cpuDebug & 1) std::cout << getName() << " integral returning " << ret << std::endl; 
  return ret; 
}

__host__ double AddPdf::sumOfNll (int numVars) const {
  static thrust::plus<double> cudaPlus;
  thrust::constant_iterator<int> eventSize(numVars); 
  thrust::constant_iterator<fptype*> arrayAddress(dev_event_array); 
  double dummy = 0;
  thrust::counting_iterator<int> eventIndex(0); 
  double ret = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)), 
					thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
					*logger, dummy, cudaPlus); 

  if (extended) {
    fptype expEvents = 0; 
    //std::cout << "Weights:"; 
    for (unsigned int i = 0; i < components.size(); ++i) {
      expEvents += host_params[host_indices[parameters + 3*(i+1)]]; 
      //std::cout << " " << host_params[host_indices[parameters + 3*(i+1)]]; 
    }
    // Log-likelihood of numEvents with expectation of exp is (-exp + numEvents*ln(exp) - ln(numEvents!)). 
    // The last is constant, so we drop it; and then multiply by minus one to get the negative log-likelihood. 
    ret += (expEvents - numEvents*log(expEvents)); 
    //std::cout << " " << expEvents << " " << numEvents << " " << (expEvents - numEvents*log(expEvents)) << std::endl; 
  }

  return ret; 
}
