#pragma once

#include <goofit/BinnedDataSet.h>
#include <goofit/PDFs/GooPdf.h>
#include <thrust/device_vector.h>

namespace GooFit {

/**
An interpolating histogram; in one dimension:

\f[
    P(x) = \frac{f(x, b(x))H[b(x)] + f(x, 1 + b(x))H[b(x) + 1]}{f(x, b(x)) + f(x, 1 + b(x))}
\f]

where \f$H\f$ is a histogram, \f$H[n]\f$ is the content of its bin with
index \f$n\f$, \f$b(x)\f$ is a function that returns the bin number that \f$x\f$
falls into, and \f$f(x, n)\f$ is the distance between \f$x\f$ and the center
of bin \f$n\f$. In other words, it does linear interpolation between
bins.

However, there are two complicating factors. First, the
histogram may have up to ten (\ref footnote21 "1") dimensions. Second, the dimensions
may be either observables or fit parameters. So, for example,
suppose we want to fit for the width \f$\sigma\f$ of a Gaussian
distribution, without using the potato of PDFs. We can do this by
making a two-dimensional histogram: The \f$x\f$ dimension is the
observable, the \f$y\f$ is \f$\sigma\f$. Fill the histogram with the value
of the Gaussian (\ref footnote22 "2") at each \f$x\f$ given the \f$\sigma\f$ in that bin. Now
when the fit asks the PDF, "What is your value at \f$x\f$ given this
\f$\sigma\f$?", the PDF responds by interpolating linearly between four
bins - ones that were precalculated with \f$\sigma\f$ values close to
what the fit is asking about. For the Gaussian this is rather
un-necessary, but may save some time for computationally expensive
functions.

The constructor takes a `BinnedDataSet` representing the underlying
histogram, a `vector` of fit parameters, and a `vector` of
observables.

\anchor footnote21 1: On the grounds that ten dimensions should be enough for anyone!

\anchor footnote22 2: Oops, there's that potato after all. It's a contrived example.
**/

class InterHistPdf : public GooPdf {
  public:
    InterHistPdf(std::string n, BinnedDataSet *x, std::vector<Variable> params, std::vector<Observable> obses);

  private:
    thrust::device_vector<fptype> *dev_base_histogram;
    fptype totalEvents;
    int numVars;
};

} // namespace GooFit
