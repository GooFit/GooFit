#include "goofit/PDFs/ScaledGaussianPdf.h"
//#include <limits>

EXEC_TARGET fptype device_ScaledGaussian (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[0]; 
  fptype mean = p[indices[1]] + p[indices[3]];
  fptype sigma = p[indices[2]] * (1 + p[indices[4]]);
  fptype ret = EXP(-0.5*(x-mean)*(x-mean)/(sigma*sigma));

#ifdef CUDAPRINT
  //if ((0 == THREADIDX) && (0 == BLOCKIDX) && (callnumber < 10)) 
    //cuPrintf("device_ScaledGaussian %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], mean, sigma, callnumber, indices, ret); 
#endif 
  //if ((gpuDebug & 1) && (0 == callnumber) && (THREADIDX == 6) && (0 == BLOCKIDX)) printf("[%i, %i] Scaled Gaussian: %f %f %f %f\n", BLOCKIDX, THREADIDX, x, mean, sigma, ret);

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
  GET_FUNCTION_ADDR(ptr_to_ScaledGaussian);
  initialise(pindices); 
}

