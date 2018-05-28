#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A variant of `AddPdf`, in which the weights
are not fit parameters but rather observables. It otherwise works
the same way as `AddPdf`; the constructor takes `vector`s of the
weights and components, and it has extended and non-extended
variants. Note that you should not mix-and-match; the weights must
be either all observables or all fit parameters.
**/

class EventWeightedAddPdf : public GooPdf {
  public:
    EventWeightedAddPdf(std::string n, std::vector<Observable> weights, std::vector<PdfBase *> comps);
    __host__ fptype normalize() const override;
    __host__ bool hasAnalyticIntegral() const override { return false; }

  protected:
    bool extended;
};
} // namespace GooFit
