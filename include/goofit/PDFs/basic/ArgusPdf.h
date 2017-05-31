#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class ArgusPdf : public GooPdf {
  public:
    ArgusPdf(std::string n, Variable *_x, Variable *m, Variable *s, bool upper, Variable *power = nullptr);
    __host__ bool hasAnalyticIntegral() const override { return false; }
    __host__ fptype integrate(fptype lo, fptype hi) const override;

  private:
};

} // namespace GooFit
