#include "ScaledGaussianPdf.hh"
//#include <limits>

__device__ fptype device_ScaledGaussian (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[0]; 
  fptype mean = p[indices[1]] + p[indices[3]];
  fptype sigma = p[indices[2]] * (1 + p[indices[4]]);
  fptype ret = EXP(-0.5*(x-mean)*(x-mean)/(sigma*sigma));

#ifdef CUDAPRINT
  //if ((0 == threadIdx.x) && (0 == blockIdx.x) && (callnumber < 10)) 
    //cuPrintf("device_ScaledGaussian %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], mean, sigma, callnumber, indices, ret); 
#endif 
  //if ((gpuDebug & 1) && (0 == callnumber) && (threadIdx.x == 6) && (0 == blockIdx.x)) printf("[%i, %i] Scaled Gaussian: %f %f %f %f\n", blockIdx.x, threadIdx.x, x, mean, sigma, ret);

  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_ScaledGaussian = device_ScaledGaussian; 

__host__ ScaledGaussianPdf::ScaledGaussianPdf (std::string n, Variable* _x, Variable* mean, Variable* sigma, Variable* delta, Variable* epsilon) 
: GooPdf(_x, n) 
{
  registerParameter(mean);
  registerParameter(sigma);
  registerParameter(delta);
  registerParameter(epsilon);

  std::vector<unsigned int> pindices;
  pindices.push_back(mean->getIndex());
  pindices.push_back(sigma->getIndex());
  pindices.push_back(delta->getIndex());
  pindices.push_back(epsilon->getIndex());
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_ScaledGaussian, sizeof(void*));
  initialise(pindices); 
}

