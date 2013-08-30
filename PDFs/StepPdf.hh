#ifndef STEP_PDF_HH
#define STEP_PDF_HH

#include "EngineCore.hh" 

class StepPdf : public EngineCore {
public:
  StepPdf (std::string n, Variable* _x, Variable* x0); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
