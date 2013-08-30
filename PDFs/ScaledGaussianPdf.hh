#ifndef SCALEDGAUSSIAN_PDF_HH
#define SCALEDGAUSSIAN_PDF_HH

#include "EngineCore.hh" 

class ScaledGaussianPdf : public EngineCore {
public:
  ScaledGaussianPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* d, Variable* e); 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 



private:

};

#endif
