#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A convolution of a classical Breit-Wigner and a
Gaussian resolution:

\f[
    P(x;m,\sigma,\Gamma) = \int\limits_{-\infty}^\infty\frac{\Gamma}{(t-m)^2-\Gamma^2/4}
e^{-\frac{(t-x)^2}{2\sigma^2}}\mathrm{d}t. \f]

The actual implementation is a horrible lookup-table-interpolation;
had Lovecraft been aware of this sort of thing, he would not have
piffled about writing about mere incomprehensible horrors from the
depths of time. The constructor takes the observable \f$x\f$, mean \f$m\f$,
Gaussian resolution width \f$\sigma\f$, and Breit-Wigner width \f$\Gamma\f$.
**/

class VoigtianPdf : public GooPdf {
  public:
    VoigtianPdf(std::string n, Observable _x, Variable m, Variable s, Variable w);
};
} // namespace GooFit
