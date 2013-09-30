#include "MappedPdf.hh"

__device__ fptype device_Mapped (fptype* evt, fptype* p, unsigned int* indices) {
  // Structure : nP mapFunctionIndex mapParamIndex functionIndex1 parameterIndex1 functionIndex2 parameterIndex2 ... 

  // Find mapping between event variables and function to evaluate
  unsigned int mapFunction = indices[1];
  // This is an index into the MappedPdf's list of functions
  //int targetFunction = (int) FLOOR(0.5 + (*(reinterpret_cast<device_function_ptr>(device_function_table[mapFunction])))(evt, p, paramIndices + indices[2]));
  int targetFunction = (int) FLOOR(0.5 + callFunction(evt, mapFunction, indices[2]));
  
  targetFunction *= 2; // Because there are two pieces of information about each function
  targetFunction += 3; // Because first function information begins at index 3

  //fptype ret = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[targetFunction]])))(evt, p, paramIndices + indices[targetFunction + 1]); 
  fptype ret = callFunction(evt, indices[targetFunction], indices[targetFunction + 1]); 
  ret *= normalisationFactors[indices[targetFunction + 1]]; 
  //if (gpuDebug & 1) 
  //if ((gpuDebug & 1) && (0 == blockIdx.x) && (0 == threadIdx.x))
    //printf("[%i, %i] Mapped: %i (%f %f %f %f) %f\n", blockIdx.x, threadIdx.x, targetFunction, evt[0], evt[1], evt[2], evt[3], ret); 
  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_Mapped = device_Mapped; 

__host__ MappedPdf::MappedPdf (std::string n, GooPdf* m, vector<GooPdf*>& t)
  : GooPdf(0, n) 
{
  components.push_back(m); 
  std::vector<unsigned int> pindices;
  pindices.push_back(m->getFunctionIndex()); 
  pindices.push_back(m->getParameterIndex()); 

  std::set<int> functionIndicesUsed;
  for (vector<GooPdf*>::iterator f = t.begin(); f != t.end(); ++f) {
    components.push_back(*f); 
    pindices.push_back((*f)->getFunctionIndex()); 
    pindices.push_back((*f)->getParameterIndex()); 
    functionIndicesUsed.insert((*f)->getFunctionIndex());
  }
  if (functionIndicesUsed.size() > 1) {
    std::cout << "Warning: More than one function type given to MappedPdf "
	      << getName()
	      << " constructor. This may slow execution by causing sequential evaluations.\n"; 
  }

  getObservables(observables); 
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_Mapped, sizeof(void*));
  initialise(pindices); 
}

__host__ fptype MappedPdf::normalise () const {
  //std::cout << "Normalising MappedPdf " << getName() << std::endl; 
  fptype ret = 0;
  for (unsigned int i = 1; i < components.size(); ++i) { // No need to normalise mapping function. 
    fptype curr = components[i]->normalise(); 
    ret += curr;
  }
  host_normalisation[parameters] = 1.0; 
  return ret; 
}
