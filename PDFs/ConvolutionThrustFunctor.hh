#ifndef CONVOLVE_THRUST_FUNCTOR_HH
#define CONVOLVE_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class ConvolutionThrustFunctor : public EngineCore {
public:

  ConvolutionThrustFunctor (std::string n, Variable* _x, EngineCore* model, EngineCore* resolution); 
  ConvolutionThrustFunctor (std::string n, Variable* _x, EngineCore* model, EngineCore* resolution, unsigned int numOthers); 
  __host__ virtual fptype normalise () const;
  __host__ void setIntegrationConstants (fptype lo, fptype hi, fptype step); 
  __host__ void registerOthers (std::vector<ConvolutionThrustFunctor*> others);

private:
  EngineCore* model;
  EngineCore* resolution; 

  fptype* host_iConsts; 
  fptype* dev_iConsts; 
  thrust::device_vector<fptype>* modelWorkSpace;
  thrust::device_vector<fptype>* resolWorkSpace; 
  int workSpaceIndex; 

};


#endif
