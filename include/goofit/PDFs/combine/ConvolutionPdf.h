#pragma once

#include <goofit/PDFs/CombinePdf.h>
#include <thrust/device_vector.h>

namespace GooFit {

/**
Numerically calculates a convolution integral

\f[
\begin{align}
    P(x;f,g) &=&
    f\otimes g \\
    &=&
    \int\limits_{-\infty}^\infty f(t) g(x-t) \mathrm{d}t.
\end{align}
\f]

The constructor takes the observable \f$x\f$, model function \f$f\f$, and
resolution function \f$g\f$.

The implementation of this function is a little complicated and
relies on caching. There is a variant constructor for cases where
several convolutions may run at the same time, eg a `MappedPdf`
where all the targets are convolutions. This variant does
cooperative loading of the caches, which is a *really neat*
optimisation and ought to work a lot better than it, actually, does.
Its constructor takes the observable, model, and resolution as
before, and an integer indicating how many other convolutions are
going to be using the same cache space.
**/

class ConvolutionPdf : public CombinePdf {
  public:
    ConvolutionPdf(std::string n, Observable _x, GooPdf *model, GooPdf *resolution);
    ConvolutionPdf(std::string n, Observable _x, GooPdf *model, GooPdf *resolution, unsigned int numOthers);
    __host__ auto normalize() -> fptype override;
    __host__ void setIntegrationConstants(fptype lo, fptype hi, fptype step);
    __host__ void registerOthers(std::vector<ConvolutionPdf *> others);

  private:
    GooPdf *model;
    GooPdf *resolution;

    fptype *host_iConsts;
    fptype *dev_iConsts;
    thrust::device_vector<fptype> *modelWorkSpace;
    thrust::device_vector<fptype> *resolWorkSpace;
    int workSpaceIndex;
};

} // namespace GooFit
