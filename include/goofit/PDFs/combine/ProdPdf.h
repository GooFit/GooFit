#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class ProdPdf : public GooPdf {
  public:
    ProdPdf(std::string n, std::vector<PdfBase *> comps);
    __host__ fptype normalize() const override;
    __host__ bool hasAnalyticIntegral() const override { return false; }

  private:
    bool varOverlaps; // True if any components share an observable.
};
} // namespace GooFit
