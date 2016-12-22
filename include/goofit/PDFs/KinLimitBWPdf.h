#ifndef KINLIMITBW_PDF_HH
#define KINLIMITBW_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class KinLimitBWPdf : public GooPdf {

public:
  KinLimitBWPdf (std::string n, Variable* _x, Variable* m, Variable* s); 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 
  __host__ void setMasses (fptype bigM, fptype smallM); 

private:

};

#endif
