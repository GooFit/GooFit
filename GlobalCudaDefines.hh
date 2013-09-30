#ifndef __GLOBAL_CUDA_HH__
#define __GLOBAL_CUDA_HH__

extern int host_callnumber; 
#include <cmath> 
#include <string> 
using namespace std; 

#ifdef OMP_ON
#include "omp.h"
#define MAX_THREADS 8
#pragma omp threadprivate (host_callnumber)
#endif

#if THRUST_DEVICE_BACKEND==THRUST_DEVICE_BACKEND_OMP
// OMP target - all 'device' memory is actually on host. 
#define MEM_DEVICE
#define MEM_SHARED
#define EXEC_TARGET __host__
#define MEMCPY(target, source, amount, dummy) memcpy(target, source, amount)
#else
// CUDA target - defaults
#define MEM_DEVICE __device__
#define MEM_SHARED __shared__
#define EXEC_TARGET __device__
#define MEMCPY(target, source, amount, direction) cudaMemcpy(target, source, amount, direction)
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
