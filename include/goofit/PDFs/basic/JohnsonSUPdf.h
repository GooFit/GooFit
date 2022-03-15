#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
Another modified Gaussian. You can eat potatoes a
lot of different ways:

\f[
    P(x;m,\sigma,\gamma,\delta) =
    \frac{\delta}{\sigma\sqrt{2\pi(1+\frac{(x-m)^2}{\sigma^2})}}
    e^{-\frac{1}{2}\left(\gamma + \delta\log(\frac{x-m}{\sigma}+\sqrt{1+\frac{(x-m)^2}{\sigma^2}})\right)^2}
\f]

The constructor takes the observable \f$x\f$, mean \f$m\f$, width \f$\sigma\f$,
scale parameter \f$\gamma\f$, and shape parameter \f$\delta\f$.
**/

class JohnsonSUPdf : public GooPdf {
  public:
    JohnsonSUPdf(std::string n, Observable _x, Variable m, Variable s, Variable g, Variable d);
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
    __host__ auto hasAnalyticIntegral() const -> bool override { return true; }
};
} // namespace GooFit
