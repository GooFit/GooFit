#ifndef __GLOBAL_CUDA_HH__
#define __GLOBAL_CUDA_HH__

#include <thrust/functional.h> // Needed for Thrust constants
#include <cmath>
#include <string>
using namespace std;
extern int host_callnumber;

//  Non-cuda defines
#if THRUST_DEVICE_SYSTEM!=THRUST_DEVICE_SYSTEM_CUDA
// OMP target - all 'device' memory is actually on host.
#define __align__(n)
#define __device__
#define __shared__
#define __constant__
#define __host__
#define DEVICE_VECTOR thrust::host_vector
// Use char* here because I need +1 to mean "offset by one byte", not "by one sizeof(whatever)".
// Can't use void* because then the compiler doesn't know how to do pointer arithmetic.
// This will fail if sizeof(char) is more than 1. But that should never happen, right?
#define MEMCPY(target, source, count, direction) memcpy((char*) target, source, count)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) memcpy(((char*) target)+offset, source, count)
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) memcpy((char*) target, ((char*) source)+offset, count)
#define GET_FUNCTION_ADDR(fname) host_fcn_ptr = (void*) fname
#define cudaDeviceSynchronize dummySynch
#define BLOCKIDX (1)
inline void cudaDeviceSynchronize() {}
#define CONST_PI M_PI
// Create my own error type to avoid __host__ redefinition
// conflict in Thrust from including driver_types.h
enum gooError {gooSuccess = 0, gooErrorMemoryAllocation};
#define RO_CACHE(x) x
#endif

#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_TBB

#include <omp.h>
#define THREADIDX (omp_get_thread_num())
#define BLOCKDIM (omp_get_num_threads())
#define THREAD_SYNCH _Pragma("omp barrier") // valid in C99 and C++11, but probably not C++93

#elif THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_CPP

#define THREADIDX (1)
#define BLOCKDIM (1)
#define THREAD_SYNCH 

#elif THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_CUDA

// CUDA target - defaults
#define THREAD_SYNCH __syncthreads();
#define DEVICE_VECTOR thrust::device_vector
#define MEMCPY(target, source, count, direction) cudaMemcpy(target, source, count, direction)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) cudaMemcpyToSymbol(target, source, count, offset, direction)
#ifdef TARGET_SM35
#define RO_CACHE(x) __ldg(&x)
#else
#define RO_CACHE(x) x
#endif
#define GET_FUNCTION_ADDR(fname) cudaMemcpyFromSymbol((void**) &host_fcn_ptr, fname, sizeof(void*))
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) cudaMemcpyFromSymbol(target, source, count, offset, direction)
// For CUDA case, just use existing errors, renamed
#include <driver_types.h>      // Needed for cudaError_t
enum gooError {gooSuccess = cudaSuccess,
               gooErrorMemoryAllocation = cudaErrorMemoryAllocation
              };
#define THREADIDX (threadIdx.x)
#define BLOCKDIM (blockDim.x)
#define BLOCKIDX (blockIdx.x)
#define CONST_PI CUDART_PI

#else

#define THREADIDX (1)
#define BLOCKDIM (1)
#endif

gooError gooMalloc(void** target, size_t bytes);
gooError gooFree(void* ptr);

#define DOUBLES 1

void abortWithCudaPrintFlush(std::string file, int line);

#ifdef DOUBLES
#define root2 1.4142135623730951
#define invRootPi 0.5641895835477563

typedef double fptype;
// Double math functions
#define ATAN2 atan2
#define COS cos
#define COSH cosh
#define SINH sinh
#define ERF erf
#define ERFC erfc
#define EXP exp
#define FABS fabs
#define FMOD fmod
#define LOG log
#define MODF modf
#define SIN sin
#define SQRT sqrt
#ifdef TARGET_SM35
#define RSQRT rsqrt
#else
#define RSQRT 1.0/SQRT
#endif
#define FLOOR floor
#define POW pow
#else
typedef float fptype;

#define root2 1.4142135623730951f
#define invRootPi 0.5641895835477563f


// Float math functions
#define ATAN2 atan2f
#define COS cosf
#define COSH coshf
#define SINH sinhf
#define ERF erff
#define ERFC erfcf
#define EXP expf
#define FABS fabsf
#define FMOD fmodf
#define LOG logf
#define MODF modff
#define SIN sinf
#define SQRT sqrtf
#define FLOOR floorf
#define POW powf
#endif


#endif
