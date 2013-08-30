#ifndef BIFURGAUSSIAN_PDF_HH
#define BIFURGAUSSIAN_PDF_HH

#include "EngineCore.hh"

class BifurGaussPdf : public EngineCore {
  public:
    BifurGaussPdf (std::string n, Variable *_x, Variable* m, Variable* sL, Variable* sR);
    __host__ fptype integrate(fptype lo, fptype hi) const;
    //__host__ virtual bool hasAnalyticIntegral () const {return true;}

  private:
 
};

#endif
