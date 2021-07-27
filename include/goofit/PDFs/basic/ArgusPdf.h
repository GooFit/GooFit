#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
Implements a threshold function

\f[
    P(x;m_0,a,p) =
    \left\{
        \begin{matrix}
        0 & x \le m_0 \\
        x\left(\frac{x^2-m_0^2}{m_0^2}\right)^p e^{a\frac{x^2-m_0^2}{m_0^2}} & x > m_0 \\
        \end{matrix}
    \right.
\f]

where the power \f$p\f$ is, by default, fixed at
0.5. The constructor takes `Variable`s representing \f$x\f$, \f$m_0\f$, and
\f$a\f$, followed by a boolean indicating whether the threshold is an
upper or lower bound. The equation above shows the PDF for a lower
bound; for upper bounds, \f$x^2-m_0^2\f$ becomes instead \f$m_0^2-x^2\f$,
and the value is zero above rather than below \f$m_0\f$. The constructor
also takes an optional `Variable` representing the power \f$p\f$; if not
given, a default parameter with value 0.5 is created.
**/

class ArgusPdf : public GooPdf {
  public:
    ArgusPdf(std::string n, Observable _x, Variable m, Variable s, bool upper);
    ArgusPdf(std::string n, Observable _x, Variable m, Variable s, bool upper, Variable power);
    __host__ auto hasAnalyticIntegral() const -> bool override { return false; }
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
};

} // namespace GooFit
