#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
An exponential decay convolved with a Gaussian
resolution:

\f[
    P(t;m,\sigma,\tau) = e^{-t/\tau} \otimes e^{-\frac{(t-m)^2}{2\sigma^2}} \\
    = (\tau/2)e^{(\tau/2)(2m+\tau\sigma^2-2t}\mathrm{erfc}\left(\frac{m+\tau\sigma^2-t}{\sigma\sqrt{2}}\right)
\f]
where \f$\mathrm{erfc}\f$ is the complementary error function. The
constructor takes the observed time \f$t\f$, mean \f$m\f$ and width \f$\sigma\f$
of the resolution, and lifetime \f$\tau\f$. Note that the original decay
function is zero for \f$t<0\f$.
**/

class ExpGausPdf : public GooPdf {
  public:
    ExpGausPdf(std::string n, Observable _x, Variable m, Variable s, Variable t);
};
} // namespace GooFit
