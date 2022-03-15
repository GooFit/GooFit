#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
Also known as the Heaviside function. Zero up to a point,
then 1 after that point:
\f[
    P(x;x_0) =
    \left\{
        \begin{matrix}
            0 & x \le x_0 \\
            1 & x > x_0
        \end{matrix}
    \right.
\f]

The constructor takes the observable \f$x\f$ and
threshold \f$x_0\f$.
**/

class StepPdf : public GooPdf {
  public:
    StepPdf(std::string n, Observable _x, Variable x0, int Up = 1);
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
    __host__ auto hasAnalyticIntegral() const -> bool override { return true; }
};
} // namespace GooFit
