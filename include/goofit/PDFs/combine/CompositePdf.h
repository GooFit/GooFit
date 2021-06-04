#pragma once

#include <goofit/PDFs/CombinePdf.h>

namespace GooFit {

/**
A chained function,

\f[
    P(x) = h(g(x)).
\f]

The constructor takes the kernel
function \f$g\f$ and the shell function \f$h\f$. Note that only
one-dimensional composites are supported - \f$h\f$ cannot take more than
one argument. The core function \f$g\f$ can take any number.

> Composites of arbitrary functions, ie \f$f(x) = h(g(x))\f$
> for any \f$h\f$ and \f$g\f$.
> In principle we should allow multi-
> dimensional compositing, eg \f$f(x, y) = i(g(x, y), h(x, y))\f$.
> Not implemented yet.
**/

class CompositePdf : public CombinePdf {
  public:
    CompositePdf(std::string n, PdfBase *core, PdfBase *shell); // Where 'core' corresponds to 'g' and 'shell' to 'h'.
    __host__ auto normalize() -> fptype override;
};
} // namespace GooFit
