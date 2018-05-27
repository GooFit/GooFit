#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A shape with a long right-hand tail - so long, in fact,
that its moments are not defined. If the most probable value (note
that this is not a mean) and the width are taken as 0 and 1, the PDF
is

\f[
P(x) = \frac{1}{\pi}\int_0^\infty e^{-t\log t - xt}\sin(t\pi)\mathrm{d}t
\f]

but the GooFit implementation is a lookup table stolen from CERNLIB.
The constructor takes the observable \f$x\f$, most probable value \f$\mu\f$
(which shifts the above expression) and the width \f$\sigma\f$ (which
scales it).
**/

class LandauPdf : public GooPdf {
  public:
    LandauPdf(std::string n, Observable _x, Variable mpv, Variable sigma);
};
} // namespace GooFit
