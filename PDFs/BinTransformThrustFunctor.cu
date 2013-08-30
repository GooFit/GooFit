#include "BinTransformThrustFunctor.hh"

__device__ fptype device_BinTransform (fptype* evt, fptype* p, unsigned int* indices) {
  // Index structure: nP lim1 bin1 lim2 bin2 ... nO o1 o2 
  int numObservables = indices[1 + indices[0]];
  int ret = 0;
  int previousSize = 1; 
  //printf("[%i, %i] Bin Transform: %i %i %f %f\n", threadIdx.x, blockIdx.x, numObservables, previousSize, evt[0], evt[1]); 
  for (int i = 0; i < numObservables; ++i) {
    fptype obsValue   = evt[indices[2 + indices[0] + i]];
    fptype lowerLimit = functorConstants[indices[i*3+1]];
    fptype binSize    = functorConstants[indices[i*3+2]];
    int numBins       = indices[i*3+3]; 

    int localBin = (int) FLOOR((obsValue - lowerLimit) / binSize);
    ret += localBin * previousSize;
    previousSize *= numBins; 
  }

  return fptype(ret); 
}

__device__ device_function_ptr ptr_to_BinTransform = device_BinTransform; 

// Notice that bin sizes and limits can be different, for this purpose, than what's implied by the Variable members. 
__host__ BinTransformThrustFunctor::BinTransformThrustFunctor (std::string n, vector<Variable*> obses, vector<fptype> limits, vector<fptype> binSizes, vector<int> numBins) 
  : ThrustPdfFunctor(0, n) 
{

  cIndex = registerConstants(2*obses.size());
  fptype* host_constants = new fptype[2*obses.size()]; 
  std::vector<unsigned int> pindices;
  for (unsigned int i = 0; i < obses.size(); ++i) {
    registerObservable(obses[i]); 
    pindices.push_back(cIndex + 2*i); 
    pindices.push_back(cIndex + 2*i + 1); 
    pindices.push_back(numBins[i]); 

    host_constants[2*i] = limits[i]; // cIndex will be accounted for by offset in memcpy
    host_constants[2*i+1] = binSizes[i]; 
  }

  cudaMemcpyToSymbol(functorConstants, host_constants, 2*obses.size()*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice); 
  delete[] host_constants; 

  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_BinTransform, sizeof(void*));
  initialise(pindices); 
}


