#ifndef __GLOBAL_CUDA_HH__
#define __GLOBAL_CUDA_HH__

#include <thrust/functional.h>
extern int host_callnumber; 
#include <cmath> 
#include <string> 
using namespace std; 

#ifdef OMP_ON
#include "omp.h"
#define MAX_THREADS 8
#pragma omp threadprivate (host_callnumber)
#endif

cudaError_t gooMalloc (void** target, size_t bytes); 
cudaError_t gooFree (void* ptr); 

#if THRUST_DEVICE_BACKEND==THRUST_DEVICE_BACKEND_OMP
// OMP target - all 'device' memory is actually on host. 
#define MEM_DEVICE
#define MEM_SHARED
#define MEM_CONSTANT 
#define EXEC_TARGET __host__
#define MEMCPY(target, source, count, dummy) memcpy(target, source, count)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) memcpy(target, source, count)
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) memcpy(target, source, count)
#define SYNCH dummySynch
void dummySynch () {}
#else
// CUDA target - defaults
#define MEM_DEVICE __device__
#define MEM_SHARED __shared__
#define MEM_CONSTANT __constant__ 
#define EXEC_TARGET __device__
#define SYNCH cudaDeviceSynchronize 
#define MEMCPY(target, source, count, direction) cudaMemcpy(target, source, count, direction) 
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) cudaMemcpyToSymbol(target, source, count, offset, direction)
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) cudaMemcpyFromSymbol(target, source, count, offset, direction)
#endif



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
