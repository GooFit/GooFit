#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class KinLimitBWPdf : public GooPdf {
  public:
    KinLimitBWPdf(std::string n, Variable *_x, Variable *m, Variable *s);
    __host__ bool hasAnalyticIntegral() const override { return false; }
    __host__ void setMasses(fptype bigM, fptype smallM);

    __host__ virtual void recursiveSetIndices ();
  private:
};
} // namespace GooFit
