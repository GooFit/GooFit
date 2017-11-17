#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ArgusPdf : public GooPdf {
  public:
    ArgusPdf(std::string n, Observable _x, Variable m, Variable s, bool upper);
    ArgusPdf(std::string n, Observable _x, Variable m, Variable s, bool upper, Variable power);
    __host__ bool hasAnalyticIntegral() const override { return false; }
    __host__ fptype integrate(fptype lo, fptype hi) const override;

    __host__ virtual void recursiveSetIndices();

  private:
    int ArgusType;
};

} // namespace GooFit
