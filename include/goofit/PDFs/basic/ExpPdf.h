#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A plain exponential:

\f[
    P(x;\alpha, x_0) = e^{\alpha(x-x_0)}
\f]

taking the
observable \f$x\f$, exponential constant \f$\alpha\f$, and optional offset
\f$x_0\f$. If \f$x_0\f$ is not specified it defaults to zero.

A variant constructor takes, in place of \f$\alpha\f$, a `vector` of coefficients
(in the order \f$\alpha_0\f$ to \f$\alpha_n\f$) to form a polynomial in the exponent:

\f[
P(x;\alpha_0, \alpha_1, \ldots \alpha_n, x_0) = e^{\alpha_0 + \alpha_1(x-x_0) + \alpha_2(x-x_0)^2 + \ldots +
\alpha_n(x-x_0)^n} \f]

The offset \f$x_0\f$ is again optional and defaults to zero.
**/

class ExpPdf : public GooPdf {
  public:
    ExpPdf(std::string n, Observable _x, Variable alpha);
    ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights);
    ExpPdf(std::string n, Observable _x, Variable alpha, Variable offset);
    ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights, Variable offset);
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
    __host__ auto hasAnalyticIntegral() const -> bool override { return (1 == host_parameters[parametersIdx]); }
};
} // namespace GooFit
