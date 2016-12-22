#ifndef CRYSTALBALL_PDF_HH
#define CRYSTALBALL_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class CrystalBallPdf : public GooPdf {
public:
  CrystalBallPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* a, Variable* power = 0); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  //__host__ virtual bool hasAnalyticIntegral () const {return true;} 


private:

};

#endif
