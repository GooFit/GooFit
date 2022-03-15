#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
What can I say? It's a normal distribution, the
potato of PDFs. Kind of bland, but goes with anything. National
cuisines have been based on it.

\f[
    P(x;m,\sigma) = e^-\frac{(x-m)^2}{2\sigma^2}
\f]

The
constructor takes the observable \f$x\f$, mean \f$m\f$, and width \f$\sigma\f$.
**/

class GaussianPdf : public GooPdf {
  public:
    GaussianPdf(std::string n, Observable _x, Variable m, Variable s);
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
    __host__ auto hasAnalyticIntegral() const -> bool override { return true; }
};
} // namespace GooFit
