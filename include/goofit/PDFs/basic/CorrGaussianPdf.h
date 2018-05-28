#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A correlated Gaussian - that is, a function of
two variables \f$x\f$ and \f$y\f$, each described by a Gaussian
distribution, but the width of the \f$y\f$ distribution depends on \f$x\f$:

\f[
    P(x,y;\bar x,\sigma_x,\bar y, \sigma_y, k) =
    e^{-\frac{(x-\bar x)^2}{2\sigma_x^2}}e^{-\frac{(y-\bar y)^2}{2(1 + k(\frac{x-\bar x}{\sigma_x})^2)\sigma_y^2}}
\f]

In other words, the effective \f$\sigma_y\f$ grows quadratically in the
normalized distance from the mean of \f$x\f$, with the quadratic term
having coefficient \f$k\f$. The constructor takes observables \f$x\f$ and
\f$y\f$, means and widths \f$\bar x\f$, \f$\sigma_x\f$, \f$\bar y\f$ and \f$\sigma_y\f$,
and coefficient \f$k\f$. Notice that if \f$k\f$ is zero, the function
reduces to a product of two Gaussians,

\f[
P(x,y;\bar x,\sigma_x,\bar y, \sigma_y) = G(x;\bar x, \sigma_x)G(y;\bar y, \sigma_y).
\f]
**/

class CorrGaussianPdf : public GooPdf {
  public:
    CorrGaussianPdf(std::string n,
                    Observable _x,
                    Observable _y,
                    Variable mean1,
                    Variable sigma1,
                    Variable mean2,
                    Variable sigma2,
                    Variable correlation);
};

} // namespace GooFit
