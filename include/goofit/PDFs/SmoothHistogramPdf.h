#pragma once

#include "goofit/PDFs/GooPdf.h"
#include "goofit/BinnedDataSet.h"

class SmoothHistogramPdf : public GooPdf {
public:
    SmoothHistogramPdf(std::string n, BinnedDataSet* x, Variable* smoothing);
    __host__ virtual fptype normalise() const;
    __host__ void extractHistogram(thrust::host_vector<fptype>& host_hist) {
        host_hist = *dev_base_histogram;
    }
    __host__ void copyHistogramToDevice(thrust::host_vector<fptype>& host_histogram);

private:
    thrust::device_vector<fptype>* dev_base_histogram;
    thrust::device_vector<fptype>* dev_smoothed_histogram;
    fptype totalEvents;
    fptype* host_constants;

    static unsigned int totalHistograms;
};
