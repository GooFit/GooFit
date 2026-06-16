#pragma once

#include <thrust/version.h>

// CCCL 3.0 (shipped with CUDA 13) removed thrust::unary_function and
// thrust::binary_function. GooFit and MCBooster functors still inherit them,
// but only for the (now unused) argument_type / result_type typedefs. Restore
// minimal definitions so the existing functors keep compiling.
//
// This header is force-included (-include) for CUDA builds so it is seen
// before any functor definition, regardless of include order.
#if THRUST_VERSION >= 300000
namespace thrust {
template <typename Argument, typename Result>
struct unary_function {
    using argument_type = Argument;
    using result_type   = Result;
};
template <typename Argument1, typename Argument2, typename Result>
struct binary_function {
    using first_argument_type  = Argument1;
    using second_argument_type = Argument2;
    using result_type          = Result;
};
} // namespace thrust
#endif
