#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class ScaledGaussianPdf : public GooPdf {
  public:
    ScaledGaussianPdf(std::string n, Variable *_x, Variable *m, Variable *s, Variable *d, Variable *e);
    __host__ bool hasAnalyticIntegral() const override { return false; }

  private:
};
} // namespace GooFit
