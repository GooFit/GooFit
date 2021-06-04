#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

class SpecialComplexSum : public thrust::binary_function<ThreeComplex, ThreeComplex, ThreeComplex> {
  public:
    __host__ __device__ auto operator()(ThreeComplex one, ThreeComplex two) -> ThreeComplex {
        return {thrust::get<0>(one) + thrust::get<0>(two),
                thrust::get<1>(one) + thrust::get<1>(two),
                thrust::get<2>(one) + thrust::get<2>(two),
                thrust::get<3>(one) + thrust::get<3>(two),
                thrust::get<4>(one) + thrust::get<4>(two),
                thrust::get<5>(one) + thrust::get<5>(two)};
    }
};

} // namespace GooFit
