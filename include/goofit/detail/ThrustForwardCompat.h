#pragma once

#include <thrust/detail/config.h> // THRUST_DEVICE_SYSTEM, THRUST_VERSION
// Modern CCCL (3.x) no longer pulls thrust::tuple / thrust::get in transitively
// through headers like <thrust/functional.h>. MCBooster (unmaintained) and some
// GooFit headers rely on that; include it once, up front, for every unit.
#include <thrust/tuple.h>

// Older Thrust defined __host__/__device__ (to nothing) for non-CUDA backends;
// modern CCCL (3.x) no longer does. On CUDA builds the toolkit still provides
// them (even in host translation units), so only define them when the device
// system is not CUDA. This header is force-included (-include) into every
// translation unit, so the macros are seen before any MCBooster/GooFit header,
// regardless of include order.
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#endif

// CCCL 3.0 (shipped with CUDA 13) removed thrust::unary_function and
// thrust::binary_function. GooFit and MCBooster functors still inherit them,
// but only for the (now unused) argument_type / result_type typedefs. Restore
// minimal definitions so the existing functors keep compiling.
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
// thrust::identity was also dropped; a couple of examples still use it.
template <typename T>
struct identity {
    using argument_type = T;
    using result_type   = T;
    __host__ __device__ constexpr auto operator()(const T &x) const -> T { return x; }
};
} // namespace thrust
#endif
