#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class CrystalBallPdf : public GooPdf {
  public:
    CrystalBallPdf(std::string n, Observable _x, Variable m, Variable s, Variable a);
    CrystalBallPdf(std::string n, Observable _x, Variable m, Variable s, Variable a, Variable power);
    __host__ fptype integrate(fptype lo, fptype hi) const override;
    //__host__ virtual bool hasAnalyticIntegral () const {return true;}

  private:
};
} // namespace GooFit
