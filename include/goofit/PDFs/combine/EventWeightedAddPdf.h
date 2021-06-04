#pragma once

#include <goofit/PDFs/CombinePdf.h>

namespace GooFit {

/**
A variant of `AddPdf`, in which the weights
are not fit parameters but rather observables. It otherwise works
the same way as `AddPdf`; the constructor takes `vector`s of the
weights and components, and it has extended and non-extended
variants. Note that you should not mix-and-match; the weights must
be either all observables or all fit parameters.
**/

class EventWeightedAddPdf : public CombinePdf {
  public:
    EventWeightedAddPdf(std::string n, std::vector<Observable> weights, std::vector<PdfBase *> comps);
    __host__ auto normalize() -> fptype override;
    __host__ auto hasAnalyticIntegral() const -> bool override { return false; }

  protected:
    bool extended;
};
} // namespace GooFit
