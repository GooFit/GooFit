#pragma once

#include <thrust/detail/config.h> // __host__, __device__ defines
#include <thrust/system_error.h>  // Error types

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>

namespace GooFit {
extern int host_callnumber;
}

#ifdef _MSC_VER
#define _Pragma(x) __pragma(x)
#endif

// Allow code to work on non-CUDA systems (beyond what is provided with thrust)
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
#define __align__(n)
inline void cudaDeviceSynchronize() {}
#define __shared__
#define __constant__
#endif

// Specialty copies
#ifdef __CUDACC__
#define MEMCPY(target, source, count, direction) cudaMemcpy(target, source, count, direction)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction)                                                     \
    cudaMemcpyToSymbol(target, source, count, offset, direction)
#define GET_FUNCTION_ADDR(fname)                                                                                       \
    {                                                                                                                  \
        cudaMemcpyFromSymbol((void **)&host_fcn_ptr, fname, sizeof(void *));                                           \
        GOOFIT_DEBUG("Using function {} in {}, {}:{}", #fname, __func__, __FILE__, __LINE__);                          \
    }
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction)                                                   \
    cudaMemcpyFromSymbol(target, source, count, offset, direction)

// This automatically selects the correct CUDA arch and expands the __ldg intrinsic to work on arbitrary types
// CUDACC only
#include <generics/ldg.h>
#define RO_CACHE(x) __ldg(&x)

#else

#define MEMCPY(target, source, count, direction) memcpy((char *)target, source, count)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) memcpy(((char *)target) + offset, source, count)
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction)                                                   \
    memcpy((char *)target, ((char *)source) + offset, count)
#define GET_FUNCTION_ADDR(fname)                                                                                       \
    {                                                                                                                  \
        host_fcn_ptr = (void *)fname;                                                                                  \
        GOOFIT_DEBUG("Using function {} in {}, {}:{}", #fname, __func__, __FILE__, __LINE__);                          \
    }
#define RO_CACHE(x) x
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP || THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_TBB
#define THREADIDX (omp_get_thread_num())
#define BLOCKDIM (omp_get_num_threads())
#define BLOCKIDX (0)
#define THREAD_SYNCH _Pragma("omp barrier")

#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP
#define THREADIDX (0)
#define BLOCKDIM (1)
#define BLOCKIDX (0)
#define THREAD_SYNCH

#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#define THREADIDX (threadIdx.x)
#define BLOCKDIM (blockDim.x)
#define BLOCKIDX (blockIdx.x)
#define THREAD_SYNCH __syncthreads();
#endif

// CUDA errors (only needed for explicit memory tranfers)
// For CUDA case, just use existing errors
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <driver_types.h>
#else
enum cudaError_t { cudaSuccess, cudaErrorMemoryAllocation };
#endif

namespace GooFit {
cudaError_t gooMalloc(void **target, size_t bytes);
cudaError_t gooFree(void *ptr);

// Allow a switch to control single vs. double precision
#ifndef GOOFIT_SINGLES

using fptype = double;
#define root2 1.4142135623730951
#define invRootPi 0.5641895835477563

#else

typedef float fptype;
#define root2 1.4142135623730951f
#define invRootPi 0.5641895835477563f

#endif
} // namespace GooFit

// Often faster than pow, and works with ints on CUDA<8
#define POW2(x) ((x) * (x))
#define POW3(x) ((x) * (x) * (x))

// Add rsqrt for everyone
#if !defined(__CUDA_ARCH__) || (__CUDA_ARCH__ < 350)
template <typename T>
__host__ __device__ T rsqrt(T val) {
    return 1.0 / sqrt(val);
}
#endif

// Fix for bug in pow(double,int) for CUDA 7 and 7.5 (device problem only)
#if defined(__CUDACC__) && __CUDACC_VER_MAJOR__ < 8
template <typename T>
__host__ __device__ T pow(T x, int y) {
    return pow(x, (T)y);
}
#endif
