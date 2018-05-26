#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A non-relativistic Breit-Wigner function, sometimes called
a Cauchy function:

\f[
    P(x;m,\Gamma) = \frac{1}{2\sqrt{\pi}}\frac{\Gamma}{(x-m)^2 + \Gamma^2/4}
\f]

The constructor takes the observable \f$x\f$, mean \f$m\f$, and width
\f$\Gamma\f$.
**/

class BWPdf : public GooPdf {
  public:
    BWPdf(std::string n, Observable _x, Variable m, Variable s);
};
} // namespace GooFit
