#pragma once

#include <thrust/functional.h>

#include "goofit/PdfBase.h"

namespace GooFit {

// Notice that operators are distinguished by the order of the operands,
// and not otherwise! It's up to the user to make his tuples correctly.
class MetricTaker : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fptype> {
  public:
    MetricTaker(PdfBase *dat, void *dev_functionPtr);
    MetricTaker(int fIdx, int pIdx);

    /// Main operator: Calls the PDF to get a predicted value, then the metric
    /// to get the goodness-of-prediction number which is returned to MINUIT.
    ///
    /// Event number, dev_event_array (pass this way for nvcc reasons), event size
    __device__ fptype operator()(thrust::tuple<int, fptype *, int> t) const;

    /// Operator for binned evaluation, no metric.
    /// Used in normalisation.
    ///
    /// Event number, event size, normalisation ranges (for binned stuff, eg integration)
    __device__ fptype operator()(thrust::tuple<int, int, fptype *> t) const;

  private:
    /// Function-pointer index of processing function, eg logarithm, chi-square, other metric.
    unsigned int metricIndex;

    /// Function-pointer index of actual PDF
    unsigned int functionIdx;

    unsigned int parameters;
};

} // namespace GooFit
