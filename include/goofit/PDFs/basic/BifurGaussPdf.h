#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class BifurGaussPdf : public GooPdf {
  public:
    BifurGaussPdf(std::string n, Variable *_x, Variable *m, Variable *sL, Variable *sR);
    __host__ fptype integrate(fptype lo, fptype hi) const override;
    //__host__ virtual bool hasAnalyticIntegral () const {return true;}

  private:
};

} // namespace GooFit
