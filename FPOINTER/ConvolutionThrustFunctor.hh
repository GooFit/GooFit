#ifndef CONVOLVE_THRUST_FUNCTOR_HH
#define CONVOLVE_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class ConvolutionThrustFunctor : public ThrustPdfFunctor {
public:

  ConvolutionThrustFunctor (std::string n, Variable* _x, ThrustPdfFunctor* model, ThrustPdfFunctor* resolution); 
  ConvolutionThrustFunctor (std::string n, Variable* _x, ThrustPdfFunctor* model, ThrustPdfFunctor* resolution, unsigned int numOthers); 
  __host__ virtual fptype normalise () const;
  __host__ void setIntegrationConstants (fptype lo, fptype hi, fptype step); 
  __host__ void registerOthers (std::vector<ConvolutionThrustFunctor*> others);

private:
  ThrustPdfFunctor* model;
  ThrustPdfFunctor* resolution; 

  fptype* host_iConsts; 
  fptype* dev_iConsts; 
  thrust::device_vector<fptype>* modelWorkSpace;
  thrust::device_vector<fptype>* resolWorkSpace; 
  int workSpaceIndex; 

};


#endif
