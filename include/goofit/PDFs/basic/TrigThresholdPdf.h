#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
Intended as part of an efficiency function,
modelling a gradual fall-off near the edges of phase space:

\f[
    P(x;a,b,t) =
    \left\{
        \begin{matrix}
            1 & d > 1/2 \\
            a + (1-a) \sin(d\pi) & \mathrm{otherwise}
        \end{matrix}
    \right.
\f]

where \f$d=b(x-t)\f$ or \f$d=b(t-x)\f$ depending on
whether the function is modelling a lower or upper threshold. The
constructor takes the observable \f$x\f$ (which will be either
\f$m^2_{12}\f$ or \f$m^2_{13}\f$), threshold value \f$t\f$, trig constant \f$b\f$,
linear constant \f$a\f$, and a boolean which if true indicates an upper
threshold. A variant constructor, for modelling a threshold in the
"third" Dalitz-plot dimension \f$m^2_{23}\f$, takes both \f$m^2_{12}\f$ and
\f$m^2_{13}\f$, and an additional mass constant \f$m\f$; it then forms
\f$x = m - m^2_{12} - m^2_{13}\f$, and otherwise does the same
calculation.
**/

class TrigThresholdPdf : public GooPdf {
  public:
    TrigThresholdPdf(
        std::string n, Observable _x, Variable thresh, Variable trigConst, Variable linConst, bool upper = true);
    TrigThresholdPdf(std::string n,
                     Observable _x,
                     Observable _y,
                     Variable thresh,
                     Variable trigConst,
                     Variable linConst,
                     Variable massConstant,
                     bool upper);
};
} // namespace GooFit
