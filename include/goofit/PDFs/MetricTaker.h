#pragma once

#include <thrust/functional.h>

#include "goofit/PdfBase.h"

class MetricTaker : public thrust::unary_function<thrust::tuple<int, fptype*, int>, fptype> {
public:

    MetricTaker(PdfBase* dat, void* dev_functionPtr);
    MetricTaker(int fIdx, int pIdx);

    /// Event number, dev_event_array (pass this way for nvcc reasons), event size
    __device__ fptype operator()(thrust::tuple<int, fptype*, int> t) const;

    /// Event number, event size, normalisation ranges (for binned stuff, eg integration)
    __device__ fptype operator()(thrust::tuple<int, int, fptype*> t) const;

private:

    /// Function-pointer index of processing function, eg logarithm, chi-square, other metric.
    unsigned int metricIndex;
    
    /// Function-pointer index of actual PDF
    unsigned int functionIdx;
    
    unsigned int parameters;


};


