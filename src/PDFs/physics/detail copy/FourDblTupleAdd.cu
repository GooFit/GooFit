#include <goofit/PDFs/physics/detail/FourDblTupleAdd.h>

namespace GooFit {

__host__ __device__ auto FourDblTupleAdd::operator()(thrust::tuple<fptype, fptype, fptype, fptype> one,
                                                     thrust::tuple<fptype, fptype, fptype, fptype> two)
    -> thrust::tuple<fptype, fptype, fptype, fptype> {
    return {thrust::get<0>(one) + thrust::get<0>(two),
            thrust::get<1>(one) + thrust::get<1>(two),
            thrust::get<2>(one) + thrust::get<2>(two),
            thrust::get<3>(one) + thrust::get<3>(two)};
}

} // namespace GooFit
