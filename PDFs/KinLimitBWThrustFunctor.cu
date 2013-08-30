#include "KinLimitBWThrustFunctor.hh"

__device__ fptype getMomentum (fptype mass, fptype pimass, fptype d0mass) {
  if (mass <= 0) return 0; 
  double lambda = mass*mass - pimass*pimass - d0mass*d0mass;
  lambda *= lambda;
  lambda -= 4*pimass*pimass*d0mass*d0mass;
  if (lambda <= 0) return 0; 
  return SQRT(0.5*lambda/mass); 
}

__device__ fptype bwFactor (fptype momentum) {
  // 2.56 = 1.6^2, comes from radius for spin-1 particle
  return 1/SQRT(1.0 + 2.56 * momentum*momentum);
}

__device__ fptype device_KinLimitBW (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]]; 
  fptype mean  = p[indices[1]];
  fptype width = p[indices[2]];
  fptype d0mass = functorConstants[indices[3]+0]; 
  fptype pimass = functorConstants[indices[3]+1]; 

  mean += d0mass;
  x += d0mass;
  
  fptype pUsingRealMass = getMomentum(mean, pimass, d0mass); 
  if (0 >= pUsingRealMass) return 0; 
  
  mean *= mean; 
  fptype pUsingX     = getMomentum(x, pimass, d0mass); 
  fptype phspfactor  = pow(pUsingX / pUsingRealMass, 3) * pow(bwFactor(pUsingX) / bwFactor(pUsingRealMass), 2); 
  fptype phspMassSq  = pow(mean - x*x, 2);
  fptype phspGammaSq = pow(width*phspfactor, 2); 

  fptype ret = (phspfactor * mean*width*width)/(phspMassSq + mean*phspGammaSq); 
#ifdef CUDAPRINT
  /*
  if (((0 == threadIdx.x) && (0 == blockIdx.x) && (callnumber < 10)) || (isnan(ret))) 
      cuPrintf("KinLimitBW %f %f %f %f %f %f %f %f %f %f\n",
	       p[indices[1]], 
	       width,
	       x - d0mass,
	       pUsingX, 
	       pUsingRealMass,
	       bwFactor(pUsingRealMass), 
	       phspfactor,
	       phspMassSq,
	       phspGammaSq,
	       ret);
  */								     
#endif

  //  if (gpuDebug & 1) printf("[%i, %i] KinLimitBW: %f %f %f %f %f\n", blockIdx.x, threadIdx.x, x, mean, width, d0mass, pimass, ret);
  return ret; 
}

__device__ device_function_ptr ptr_to_KinLimitBW = device_KinLimitBW; 

__host__ KinLimitBWThrustFunctor::KinLimitBWThrustFunctor (std::string n, Variable* _x, Variable* mean, Variable* width) 
: EngineCore(_x, n) 
{
  registerParameter(mean);
  registerParameter(width);

  std::vector<unsigned int> pindices;
  pindices.push_back(mean->getIndex());
  pindices.push_back(width->getIndex());
  pindices.push_back(registerConstants(2));
  setMasses(1.8645, 0.13957); 
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_KinLimitBW, sizeof(void*));
  initialise(pindices);
}

__host__ void KinLimitBWThrustFunctor::setMasses (fptype bigM, fptype smallM) {
  fptype constants[2];
  constants[0] = bigM;
  constants[1] = smallM;
  cudaMemcpyToSymbol(functorConstants, constants, 2*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice); 
}
