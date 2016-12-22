#ifndef JOHNSONSU_PDF_HH
#define JOHNSONSU_PDF_HH

#include "GooPdf.hh" 

class JohnsonSUPdf : public GooPdf {
public:
  JohnsonSUPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* g, Variable* d); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
