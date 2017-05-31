#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class JohnsonSUPdf : public GooPdf {
  public:
    JohnsonSUPdf(std::string n, Variable *_x, Variable *m, Variable *s, Variable *g, Variable *d);
    __host__ fptype integrate(fptype lo, fptype hi) const override;
    __host__ bool hasAnalyticIntegral() const override { return true; }

  private:
};
} // namespace GooFit
