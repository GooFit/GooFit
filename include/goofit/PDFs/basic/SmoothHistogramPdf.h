#pragma once

#include "goofit/BinnedDataSet.h"
#include "goofit/PDFs/GooPdf.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace GooFit {

class SmoothHistogramPdf : public GooPdf {
  public:
    SmoothHistogramPdf(std::string n, BinnedDataSet *x, Variable *smoothing);
    __host__ fptype normalize() const override;
    __host__ void extractHistogram(thrust::host_vector<fptype> &host_hist) { host_hist = *dev_base_histogram; }
    __host__ void copyHistogramToDevice(thrust::host_vector<fptype> &host_histogram);

  private:
    thrust::device_vector<fptype> *dev_base_histogram;
    thrust::device_vector<fptype> *dev_smoothed_histogram;
    fptype totalEvents;
    fptype *host_constants;

    static unsigned int totalHistograms;
};
} // namespace GooFit
