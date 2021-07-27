#pragma once

#include <goofit/BinnedDataSet.h>
#include <goofit/PDFs/GooPdf.h>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace GooFit {

/**
Another histogram, but this one does smoothing
in place of interpolation. That is, suppose the event falls in bin
\f$N\f$ of a one-dimensional histogram; then the returned value is a
weighted average of bins \f$N-1\f$, \f$N\f$, and \f$N+1\f$. For multidimensional
cases the weighted average is over all the neighbouring bins,
including diagonals:

\f[
    P(\vec x;s;H) = \frac{H(\mathrm{bin}(\vec x)) + s\sum\limits_{i=\mathrm{neighbours}}\delta{i}H(i)}{1 +
s\sum\limits_{i=\mathrm{neighbours}}\delta{i}} \f]

where \f$\delta_i\f$ is zero for bins that fall outside the histogram
limits, and one otherwise. The constructor takes the underlying
histogram \f$H\f$ (which also defines the event vector \f$\vec x\f$) and the
smoothing factor \f$s\f$; notice that if \f$s\f$ is zero, the PDF reduces to
a simple histogram lookup. The `BinnedDataSet` representing \f$H\f$ may
be empty; in that case the lookup table should be set later using
the `copyHistogramToDevice` method.
**/

class SmoothHistogramPdf : public GooPdf {
  public:
    SmoothHistogramPdf(std::string n, BinnedDataSet *x, Variable smoothing);
    __host__ auto normalize() -> fptype override;
    __host__ void extractHistogram(thrust::host_vector<fptype> &host_hist) { host_hist = *dev_base_histogram; }
    __host__ void copyHistogramToDevice(thrust::host_vector<fptype> &host_histogram);

  private:
    thrust::device_vector<fptype> *dev_base_histogram;
    thrust::device_vector<fptype> *dev_smoothed_histogram;
    fptype totalEvents;

    static unsigned int totalHistograms;
};
} // namespace GooFit
