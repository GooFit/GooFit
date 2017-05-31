#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class AddPdf : public GooPdf {
  public:
    AddPdf(std::string n, std::vector<Variable *> weights, std::vector<PdfBase *> comps);
    AddPdf(std::string n, Variable *frac1, PdfBase *func1, PdfBase *func2);
    __host__ fptype normalize() const override;
    __host__ bool hasAnalyticIntegral() const override { return false; }

  protected:
    __host__ double sumOfNll(int numVars) const override;

  private:
    bool extended;
};

} // namespace GooFit
