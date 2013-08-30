#ifndef INTERHIST_PDF_HH
#define INTERHIST_PDF_HH

#include "EngineCore.hh" 
#include "BinnedDataSet.hh" 

class InterHistPdf : public EngineCore {
public:
  InterHistPdf (std::string n, 
			  BinnedDataSet* x, 
			  std::vector<Variable*> params, 
			  std::vector<Variable*> obses); 
  //__host__ virtual fptype normalise () const;

private:
  thrust::device_vector<fptype>* dev_base_histogram; 
  fptype totalEvents; 
  fptype* host_constants;
  int numVars; 
};

#endif
