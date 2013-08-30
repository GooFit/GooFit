#include "ExpGausPdf.hh"

__device__ fptype device_ExpGaus (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype mean  = p[indices[1]];
  fptype sigma = p[indices[2]];
  fptype alpha = p[indices[3]];

  fptype ret = 0.5*alpha; 
  fptype exparg = ret * (2*mean + alpha*sigma*sigma - 2*x);
  fptype erfarg = (mean + alpha*sigma*sigma - x) / (sigma * 1.4142135623);

  ret *= EXP(exparg); 
  ret *= ERFC(erfarg); 

  return ret; 
}

__device__ device_function_ptr ptr_to_ExpGaus = device_ExpGaus; 

ExpGausPdf::ExpGausPdf (std::string n, Variable* _x, Variable* mean, Variable* sigma, Variable* tau) 
  : EngineCore(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));
  pindices.push_back(registerParameter(tau));
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_ExpGaus, sizeof(void*));
  initialise(pindices); 
}


