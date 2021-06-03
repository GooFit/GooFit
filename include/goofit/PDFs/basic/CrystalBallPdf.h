#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A Gaussian with a power-law tail on one side:

\f[
\begin{align}
    P(x;m,\sigma,\alpha,p) &=& \left\{ \begin{matrix}
    e^{-\frac{(x-m)^2}{2\sigma^2}} & \mathrm{sg}(\alpha)\frac{x - m}{\sigma} \le \mathrm{sg}(\alpha)\alpha \\
    e^{-\alpha^2/2}\left(\frac{p/\alpha}{p/\alpha - \alpha + \frac{x-m}{\sigma}}\right)^p
    & \mathrm{otherwise } (\alpha\ne 0). \\
    \end{matrix}
    \right.
\end{align}
\f]

The constructor takes the observable \f$x\f$,
the mean \f$m\f$, width \f$\sigma\f$, cutoff \f$\alpha\f$, and power \f$p\f$. Note
that if \f$\alpha\f$ is negative, the power-law tail is on the right; if
positive, on the left. For \f$\alpha=0\f$, the function reduces to a
simple Gaussian in order to avoid \f$p/\alpha\f$ blowing up.
**/

class CrystalBallPdf : public GooPdf {
  public:
    CrystalBallPdf(std::string n, Observable _x, Variable m, Variable s, Variable a);
    CrystalBallPdf(std::string n, Observable _x, Variable m, Variable s, Variable a, Variable power);
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
    //__host__ virtual bool hasAnalyticIntegral () const {return true;}

  private:
};
} // namespace GooFit
