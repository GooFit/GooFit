#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

const std::string ArgusHelp = R"raw(
Implements a threshold function

$$
    P(x;m_0,a,p) =
    \left\{
        \begin{matrix}
            0 & x \le m_0 \\
            x\left(\frac{x^2-m_0^2}{m_0^2}\right)^p e^{a\frac{x^2-m_0^2}{m_0^2}} & x > m_0 \\
        \end{matrix}
    \right.
$$

where the power $p$ is, by default, fixed at
0.5. The constructor takes `Variable`s representing $x$, $m_0$, and
$a$, followed by a boolean indicating whether the threshold is an
upper or lower bound. The equation above shows the PDF for a lower
bound; for upper bounds, $x^2-m_0^2$ becomes instead $m_0^2-x^2$,
and the value is zero above rather than below $m_0$. The constructor
also takes an optional `Variable` representing the power $p$; if not
given, a default parameter with value 0.5 is created.
)raw";

class ArgusPdf : public GooPdf {
  public:
    ArgusPdf(std::string n, Observable _x, Variable m, Variable s, bool upper);
    ArgusPdf(std::string n, Observable _x, Variable m, Variable s, bool upper, Variable power);
    __host__ bool hasAnalyticIntegral() const override { return false; }
    __host__ fptype integrate(fptype lo, fptype hi) const override;
};

} // namespace GooFit
