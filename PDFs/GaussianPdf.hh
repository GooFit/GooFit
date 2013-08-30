#ifndef GAUSSIAN_PDF_HH
#define GAUSSIAN_PDF_HH

#include "EngineCore.hh" 

class GaussianPdf : public EngineCore {
public:
  GaussianPdf (std::string n, Variable* _x, Variable* m, Variable* s); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
