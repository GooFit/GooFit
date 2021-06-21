#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A two-sided Gaussian, with a \f$\sigma\f$ that varies
    depending on which side of the mean you are on:

\f[
    P(x;m,\sigma_L,\sigma_R) =
    \left\{
        \begin{matrix}
            e^{-\frac{(x-m)^2}{2\sigma_L^2}} & x \le m \\
            e^{-\frac{(x-m)^2}{2\sigma_R^2}} & x > m. \\
        \end{matrix}
    \right.
\f]

The constructor takes the observable \f$x\f$,
mean \f$m\f$, and left and right sigmas \f$\sigma_{L,R}\f$.
**/

class BifurGaussPdf : public GooPdf {
  public:
    BifurGaussPdf(std::string n, Observable _x, Variable m, Variable sL, Variable sR);
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
    //__host__ virtual bool hasAnalyticIntegral () const {return true;}

  private:
};

} // namespace GooFit
