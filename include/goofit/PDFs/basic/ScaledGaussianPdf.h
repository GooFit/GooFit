#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
Another Gaussian variant. This one moves its
mean by a bias \f$b\f$ and scales its width by a scale factor
\f$\epsilon\f$:

\f[
    P(x;m,\sigma,b,\epsilon) = e^{-\frac{(x+b-m)^2}{2(\sigma(1+\epsilon))^2}}.
\f]

This has a somewhat specialised function: It allows fitting Monte
Carlo to, for example, a sum of two Gaussians, whose means and
widths are then frozen. Then real data can be fit for a common bias
and \f$\epsilon\f$.

The constructor takes the observable \f$x\f$, mean \f$m\f$, width \f$\sigma\f$,
bias \f$b\f$ and scale factor \f$\epsilon\f$.
**/

class ScaledGaussianPdf : public GooPdf {
  public:
    ScaledGaussianPdf(std::string n, Observable _x, Variable m, Variable s, Variable d, Variable e);
    __host__ auto hasAnalyticIntegral() const -> bool override { return false; }
};
} // namespace GooFit
