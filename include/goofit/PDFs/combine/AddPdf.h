#pragma once

#include <goofit/PDFs/CombinePdf.h>

namespace GooFit {

/**
A weighted sum of two or more PDFs. There are two
variants, `extended` and `unextended`. In the extended version the
weights are interpreted as numbers of events, and \f$N\f$ PDFs have \f$N\f$
weights; in the unextended version the weights are probabilities
(i.e., between 0 and 1) and \f$N\f$ PDFs have \f$N-1\f$ weights, with the
probability of the last PDF being 1 minus the sum of the weights of
the others.

\f[
\begin{align}
    P(F_1,\ldots, F_n,w_1,\ldots,w_n)
    &=&
    w_1F_1 + \ldots + w_nF_n \\
    P(F_1,\ldots, F_n,w_1,\ldots,w_{n-1})
    &=&
    w_1F_1 + \ldots + w_{n-1}F_{n-1}\\
    &&
    + (1 - w_1 - \ldots - w_{n-1})F_n.
\end{align}
\f]

The constructor
takes a `vector` of weights \f$w_i\f$ and a `vector` of components
\f$F_i\f$. If the two `vector`s are of equal length the extended version
is used; if there is one more component than weight, the unextended
version; anything else is an error. There is also a special-case
constructor taking a single weight and two components, to save
creating the `vector`s in this common case.

Note that this PDF overrides the `sumOfNll` method; if an extended
`AddPdf` is used as a top-level PDF (that is, sent to `FitManager`
for fitting), an additional term for the number of events will be
added to the NLL.

Also note that if the `AddPdf` has `setCommonNorm()` called, the normalization
changes. By default the components are normalized separately, so
that

\f[
P(x;\vec F, \vec w) =
\sum\limits_i \frac{w_iF_i(x)}{\int F_i(x) \mathrm{d}x},
\f]

but with `ForceCommonNorm` set, the integral is instead taken at the
level of the sum:

\f[
P(x;\vec F, \vec w) =
\frac{\sum\limits_i w_iF_i(x)}{\int\sum\limits_i w_iF_i(x)\mathrm{d}x}.
\f]

The difference is subtle but sometimes important.
**/

class AddPdf : public CombinePdf {
  public:
    AddPdf(std::string n, std::vector<Variable> weights, std::vector<PdfBase *> comps);
    AddPdf(std::string n, Variable frac1, PdfBase *func1, PdfBase *func2);
    __host__ auto normalize() -> fptype override;
    __host__ auto hasAnalyticIntegral() const -> bool override { return false; }

  protected:
    __host__ auto calculateNLL() -> double override;

  private:
    bool extended;
};

} // namespace GooFit
