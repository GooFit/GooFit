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
    __host__ __device__ auto operator()(thrust::tuple<fptype, fptype, fptype, fptype> one,
                                        thrust::tuple<fptype, fptype, fptype, fptype> two)
        -> thrust::tuple<fptype, fptype, fptype, fptype>;
};

} // namespace GooFit
