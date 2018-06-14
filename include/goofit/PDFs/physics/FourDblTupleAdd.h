#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <mcbooster/GTypes.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class FourDblTupleAdd : public thrust::binary_function<thrust::tuple<fptype, fptype, fptype, fptype>,
                                                       thrust::tuple<fptype, fptype, fptype, fptype>,
                                                       thrust::tuple<fptype, fptype, fptype, fptype>> {
  public:
    __host__ __device__ thrust::tuple<fptype, fptype, fptype, fptype>
    operator()(thrust::tuple<fptype, fptype, fptype, fptype> one, thrust::tuple<fptype, fptype, fptype, fptype> two) {
        return {thrust::get<0>(one) + thrust::get<0>(two),
                thrust::get<1>(one) + thrust::get<1>(two),
                thrust::get<2>(one) + thrust::get<2>(two),
                thrust::get<3>(one) + thrust::get<3>(two)};
    }
};

} // namespace GooFit
