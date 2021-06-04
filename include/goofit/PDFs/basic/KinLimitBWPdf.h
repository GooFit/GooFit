#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A relativistic Breit-Wigner function modified by a
factor accounting for limited phase space (\ref footnote31 "1"); for example, in the
decay \f$D^{*+}\to D^0\pi^+\f$, the difference between the \f$D^*\f$ and
\f$D^0\f$ masses is only slightly more than the pion mass. Consequently,
the distribution of \f$\Delta m = m(D^*) - m(D^0)\f$ is slightly
asymmetric: The left side of the peak, where the phase space narrows
rapidly, is less likely than the right side.

\f[
\begin{align}
    P(x;x_0,\Gamma,M,m) &=& \left\{ \begin{matrix}
    0 & \lambda(x_0,M,m) \le 0 \\
    \frac{S(x,x_0,M,m)x_0'\Gamma^2}{\left(x_0'-x'^2\right)^2 + x_0'\Gamma^2S^2(x,x_0,M,m)} & \mathrm{otherwise.}
    \end{matrix}
    \right.
\end{align}
\f]

Here priming indicates addition of \f$M\f$, so
that \f$x'=x+M\f$, \f$x_0'=x_0+M\f$; the phase-space function \f$S\f$ and its
supporting characters \f$\lambda\f$, \f$p\f$, and \f$b_W\f$ are given by

\f[
\begin{align}
    S(x,x_0,M,m)   &=& \left(\frac{p(x,M,m)}{p(x_0,M,m)}\right)^3\left(\frac{b_W(x,M,m)}{b_W(x_0,M,m)}\right)^2 \\
    b_W(x,M,m)     &=& \frac{1}{\sqrt{1 + r^2p^2(x,M,m)}}\\
    p(x,M,m)       &=& \sqrt{\lambda(x,M,m)/(2x)}\\
    \lambda(x,M,m) &=& \left(x'^2-(M-m)^2\right)\left(x'^2-(M+m)^2\right).
\end{align}
\f]

The radius \f$r\f$ that appears in \f$b_W\f$ (which does not stand for
Breit-Wigner, but Blatt-Weisskopf!) is hardcoded to be 1.6.

The constructor takes the observable \f$x\f$, mean \f$x_0\f$, and width
\f$\Gamma\f$. The large and small masses \f$M\f$ and \f$m\f$, which determine
the phase space, are by default 1.8645 (the \f$D^0\f$ mass) and 0.13957
(mass of a charged pion), but can be set with a call to `setMasses`.
Note that they are constants, not fit parameters.

\anchor footnote31 1: If this seems complicated, spare a thought for the hapless
ergrad who had to code the original CPU version.

**/

class KinLimitBWPdf : public GooPdf {
  public:
    KinLimitBWPdf(std::string n, Observable _x, Variable m, Variable s);
    __host__ auto hasAnalyticIntegral() const -> bool override { return false; }
    __host__ void setMasses(fptype bigM, fptype smallM);
};
} // namespace GooFit
