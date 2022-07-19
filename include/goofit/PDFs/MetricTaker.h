#pragma once

#include <thrust/functional.h>

#include <goofit/PdfBase.h>

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#define cuda_error_check(stat)                                                                                         \
    { cudaErrCheck_((stat), __FILE__, __LINE__); }

void inline cudaErrCheck_(cudaError_t stat, const char *file, int line) {
    if(stat != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(stat), file, line);
    }
}
#endif

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
    __device__ auto operator()(thrust::tuple<int, fptype *, int> t) const -> fptype;

    /// Operator for binned evaluation, no metric.
    /// Used in normalization.
    ///
    /// Event number, event size, normalization ranges (for binned stuff, eg integration)
    __device__ auto operator()(thrust::tuple<int, int, fptype *> t) const -> fptype;

    /// Update which index we need to use:
    __host__ void setFunctionIndex(const int &id) { functionIdx = id; }

  private:
    /// Function-pointer index of processing function, eg logarithm, chi-square, other metric.
    unsigned int metricIndex;

    /// Function-pointer index of actual PDF
    unsigned int functionIdx;

    unsigned int parameters;
};

} // namespace GooFit
