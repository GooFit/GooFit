#ifndef CONVOLVE_PDF_HH
#define CONVOLVE_PDF_HH

#include "EngineCore.hh" 

class ConvolutionPdf : public EngineCore {
public:

  ConvolutionPdf (std::string n, Variable* _x, EngineCore* model, EngineCore* resolution); 
  ConvolutionPdf (std::string n, Variable* _x, EngineCore* model, EngineCore* resolution, unsigned int numOthers); 
  __host__ virtual fptype normalise () const;
  __host__ void setIntegrationConstants (fptype lo, fptype hi, fptype step); 
  __host__ void registerOthers (std::vector<ConvolutionPdf*> others);

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
