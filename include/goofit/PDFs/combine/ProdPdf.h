#pragma once

#include <goofit/PDFs/CombinePdf.h>

namespace GooFit {

/**
A product of two or more PDFs:
\f[
    P(x; \vec F) = \prod\limits_i F_i(x).
\f]

The
constructor just takes a `vector` of the functions to be multiplied.

GooFit::ProdPdf does allow variable overlaps, that is, the components may
depend on the same variable, eg \f$P(x) = A(x)B(x)\f$. If this happens,
the entire GooFit::ProdPdf object will be normalized together, since in
general
\f$\int A(x)B(x) \mathrm{d}x \ne \int A(x) \mathrm{d}x \int B(x) \mathrm{d}x\f$.
However, if any of the components have the flag `ForceSeparateNorm`
set, as well as in the default case that the components depend on
separate observables, each component will be normalized
individually. Some care is indicated when using the
`ForceSeparateNorm` flag, and possibly a rethink of why there is a
product of two PDFs depending on the same variable in the first
place.
**/

class ProdPdf : public CombinePdf {
  public:
    ProdPdf(std::string n, std::vector<PdfBase *> comps);
    __host__ auto normalize() -> fptype override;
    __host__ auto hasAnalyticIntegral() const -> bool override { return false; }

  private:
    bool varOverlaps; // True if any components share an observable.
};
} // namespace GooFit
