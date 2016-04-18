#ifndef __GLOBAL_CUDA_HH__
#define __GLOBAL_CUDA_HH__

#include <thrust/functional.h> // Needed for Thrust constants
#include <cmath> 
#include <string> 
using namespace std; 
extern int host_callnumber; 

// Thrust 1.7 will make the use of THRUST_DEVICE_BACKEND an error
#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP
#if THRUST_VERSION < 100699
// Ensure backwards compatibility with Thrust 1.5
#undef THRUST_DEVICE_BACKEND
#define THRUST_DEVICE_BACKEND THRUST_DEVICE_BACKEND_OMP
#endif 
// OMP target - all 'device' memory is actually on host. 
#define ALIGN(n)
#define MEM_DEVICE
#define MEM_SHARED
#define MEM_CONSTANT 
#define EXEC_TARGET __host__
#define THREAD_SYNCH _Pragma("omp barrier") // valid in C99 and C++11, but probably not C++93
#define DEVICE_VECTOR thrust::host_vector
// Use char* here because I need +1 to mean "offset by one byte", not "by one sizeof(whatever)".
// Can't use void* because then the compiler doesn't know how to do pointer arithmetic. 
// This will fail if sizeof(char) is more than 1. But that should never happen, right? 
#define MEMCPY(target, source, count, direction) memcpy((char*) target, source, count)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) memcpy(((char*) target)+offset, source, count)
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) memcpy((char*) target, ((char*) source)+offset, count)
#define GET_FUNCTION_ADDR(fname) host_fcn_ptr = (void*) fname
#define SYNCH dummySynch
#define THREADIDX (omp_get_thread_num())
#define BLOCKDIM (omp_get_num_threads())
#define BLOCKIDX (1)
void dummySynch (); 
#define CONST_PI M_PI
// Create my own error type to avoid __host__ redefinition
// conflict in Thrust from including driver_types.h
enum gooError {gooSuccess = 0, gooErrorMemoryAllocation};
#else
// CUDA target - defaults
#define ALIGN(n) __align__(n)
#define MEM_DEVICE __device__
#define MEM_SHARED __shared__
#define MEM_CONSTANT __constant__ 
#define EXEC_TARGET __device__
#define SYNCH cudaDeviceSynchronize 
#define THREAD_SYNCH __syncthreads(); 
#define DEVICE_VECTOR thrust::device_vector
#define MEMCPY(target, source, count, direction) cudaMemcpy(target, source, count, direction) 
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) cudaMemcpyToSymbol(target, source, count, offset, direction)
#define GET_FUNCTION_ADDR(fname) cudaMemcpyFromSymbol((void**) &host_fcn_ptr, fname, sizeof(void*))
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) cudaMemcpyFromSymbol(target, source, count, offset, direction)
// For CUDA case, just use existing errors, renamed
#include <driver_types.h>      // Needed for cudaError_t
enum gooError {gooSuccess = cudaSuccess, 
	       gooErrorMemoryAllocation = cudaErrorMemoryAllocation};
#define THREADIDX (threadIdx.x)
#define BLOCKDIM (blockDim.x)
#define BLOCKIDX (blockIdx.x)
#define CONST_PI CUDART_PI
#endif

gooError gooMalloc (void** target, size_t bytes); 
gooError gooFree (void* ptr); 

#define DOUBLES 1

void abortWithCudaPrintFlush (std::string file, int line); 

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
