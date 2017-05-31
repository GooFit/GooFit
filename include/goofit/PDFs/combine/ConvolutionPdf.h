#pragma once

#include "goofit/PDFs/GooPdf.h"
#include <thrust/device_vector.h>

namespace GooFit {

class ConvolutionPdf : public GooPdf {
  public:
    ConvolutionPdf(std::string n, Variable *_x, GooPdf *model, GooPdf *resolution);
    ConvolutionPdf(std::string n, Variable *_x, GooPdf *model, GooPdf *resolution, unsigned int numOthers);
    __host__ fptype normalize() const override;
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
