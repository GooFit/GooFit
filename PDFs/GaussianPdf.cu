#include "GaussianPdf.hh"

__device__ fptype device_Gaussian (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]]; 
  fptype mean = p[indices[1]];
  fptype sigma = p[indices[2]];

  fptype ret = EXP(-0.5*(x-mean)*(x-mean)/(sigma*sigma));

  //if ((0 == threadIdx.x) && (0 == blockIdx.x)) cuPrintf("Gaussian Values %f %i %i %f %f %i\n", x, indices[1], indices[2], mean, sigma, callnumber); 
  //cuPrintf("device_Gaussian %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], mean, sigma, callnumber, indices, ret); 
  //if ((0 == threadIdx.x) && (0 == blockIdx.x))
  //printf("device_Gaussian %f %f %f %i %f\n", x, mean, sigma, callnumber, ret);     


  return ret; 
}

__device__ device_function_ptr ptr_to_Gaussian = device_Gaussian; 

__host__ GaussianPdf::GaussianPdf (std::string n, Variable* _x, Variable* mean, Variable* sigma) 
  : EngineCore(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_Gaussian, sizeof(void*));
  initialise(pindices); 
}

__host__ fptype GaussianPdf::integrate (fptype lo, fptype hi) const {
  //static const fptype root2 = sqrt(2.);
  static const fptype rootPi = sqrt(atan2(0.0,-1.0));
  //static const fptype rootPiBy2 = rootPi / root2;
  
  unsigned int* indices = host_indices+parameters; 
  //fptype xscale = root2*host_params[indices[2]];

  /*
  std::cout << "Gaussian integral: " 
	    << xscale << " "
	    << host_params[indices[1]] << " "
	    << host_params[indices[2]] << " "
	    << ERF((hi-host_params[indices[1]])/xscale) << " "
	    << ERF((lo-host_params[indices[1]])/xscale) << " "
	    << rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) -
						  ERF((lo-host_params[indices[1]])/xscale)) 
	    << std::endl; 
  */
  //return rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) - 
  //					    ERF((lo-host_params[indices[1]])/xscale));

  // Integral over all R. 
  fptype sigma = host_params[indices[2]];
  sigma *= root2*rootPi;
  return sigma; 
}

