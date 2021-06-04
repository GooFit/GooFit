#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
If the Gaussian is the potato, what is the
polynomial? Bread? Milk? Nothing exotic, at any rate. The GooFit
version does have some subtleties, to allow for polynomials over an
arbitrary number (\ref footnote41 "1") of dimensions:

\f[
    P(\vec x; \vec a, \vec x_0, N) =
    \sum\limits_{p_1+p_2+\ldots+p_n \le N} a_{p_1p_2\ldots p_n} \prod\limits_{i=1}^n (\vec x - \vec x_0)_i^{p_i}
\f]

where \f$N\f$ is the highest degree of the polynomial and \f$n\f$ is the
number of dimensions. The constructor takes a `vector` of
observables, denoted \f$\vec x\f$ above; a `vector` of coefficients,
\f$\vec a\f$, a `vector` of optional offsets \f$\vec x_0\f$ (if not
specified, these default to zero), and the maximum degree \f$N\f$. The
coefficients are in the order
\f$a_{p_0p_0\ldots p_0}, a_{p_1p_0\ldots p_0}, \ldots a_{p_Np_0\ldots p_0}, a_{p_0p_1\ldots p_0}, a_{p_1p_1\ldots p_0},
\ldots a_{p_0p_0\ldots p_N}\f$. In other words, start at the index for
the constant term, and increment the power of the leftmost
observable. Every time the sum of the powers reaches \f$N\f$, reset the
leftmost power to zero and increment the next-leftmost. When the
next-leftmost reaches \f$N\f$, reset it to zero and increment the
third-leftmost, and so on.

An example may be helpful; for two
dimensions \f$x\f$ and \f$y\f$, and a maximum power of 3, the order is
\f$a_{00}, a_{10}, a_{20}, a_{30}, a_{01}, a_{11}, a_{21}, a_{02}, a_{12}, a_{03}\f$.
This can be visualised as picking boxes out of a matrix and
discarding the ones where the powers exceed the maximum:

\f[
\begin{array}{cccc}
    9: x^0y^3 &    -      &    -      &    -      \\
    7: x^0y^2 & 8: x^1y^2 &    -      &    -      \\
    4: x^0y^1 & 5: x^1y^1 & 6: x^2y^1 &    -      \\
    0: x^0y^0 & 1: x^1y^0 & 2: x^2y^0 & 3: x^3y^0 \\
\end{array}
\f]

starting in the lower-lefthand corner and going right,
then up.

There is also a simpler version of the constructor for the case of a
polynomial with only one dimension; it takes the observable, a
`vector` of coefficients, an optional offset, and the lowest (not
highest) degree of the polynomial; the latter two both default to
zero. In this case the order of the coefficients is from lowest to
highest power.

\anchor footnote41 1: Although being honest, just supporting the special cases of one and two would likely have
sufficed.
**/

class PolynomialPdf : public GooPdf {
  public:
    PolynomialPdf(std::string n, Observable _x, std::vector<Variable> weights, unsigned int lowestDegree = 0);
    PolynomialPdf(
        std::string n, Observable _x, std::vector<Variable> weights, Variable x0, unsigned int lowestDegree = 0);
    PolynomialPdf(std::string n,
                  std::vector<Observable> obses,
                  std::vector<Variable> coeffs,
                  std::vector<Variable> offsets,
                  unsigned int maxDegree = 0);
    __host__ auto integrate(fptype lo, fptype hi) const -> fptype override;
    __host__ auto getCoefficient(int coef) -> fptype;

  private:
    std::unique_ptr<Variable> center;
};
} // namespace GooFit
