#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A custom shape with a long tail:

\f[
    P(x;m,\sigma,t) =
    e^{-\frac{1}{2}\left(\log^2(1+t\frac{x-m}{\sigma}\frac{\sinh(t\sqrt{\log(4)})}{\sqrt{\log(4)}})/t + t^2\right)}
\f]

The constructor takes the observable \f$x\f$, mean \f$m\f$, width \f$\sigma\f$,
and tail factor \f$t\f$. If \f$t\f$ is less than \f$10^{-7}\f$, the function
returns a simple Gaussian, which probably indicates that it
approximates a Gaussian for small tail parameters, but I'd hate to
have to show such a thing.
**/

class NovosibirskPdf : public GooPdf {
  public:
    NovosibirskPdf(std::string n, Observable _x, Variable m, Variable s, Variable t);
};
} // namespace GooFit
