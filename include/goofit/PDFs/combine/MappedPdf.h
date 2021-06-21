#pragma once

#include <goofit/PDFs/CombinePdf.h>

namespace GooFit {

/**
A function having the form:

\f[
    F(x) =
    \left\{
        \begin{matrix}
            F_1(x)   & x_0 \le x \le x_1   \\
            F_2(x)   & x_1 < x \le x_2     \\
            (\ldots) & (\ldots)            \\
            F_n(x)   & x_{n-1} < x \le x_n \\
        \end{matrix}
    \right.
\f]

The constructor takes a *mapping function*
\f$m\f$, which returns an index; and a `vector` of evaluation functions
\f$\vec F\f$, so that if \f$m\f$ is zero, the PDF returns \f$F_0\f$, and so on.
Notice that \f$m\f$ does not strictly need to return an integer - in
fact the constraints of GooFit force it to return a floating-point
number - since GooFit::MappedPdf will round the result to the nearest
whole number. The canonical example of a mapping function is
GooFit::BinTransformPdf.
**/

class MappedPdf : public CombinePdf {
  public:
    MappedPdf(std::string n, GooPdf *m, std::vector<GooPdf *> &t);
    // Map function m must be custom written to correspond to order of function list t.
    __host__ auto normalize() -> fptype override;
};
} // namespace GooFit
