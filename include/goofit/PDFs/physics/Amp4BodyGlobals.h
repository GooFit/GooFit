#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

extern __constant__ unsigned int AmpIndices[500];
extern __device__ fpcomplex *cResSF[10];
extern __device__ fpcomplex *cResSF_TD[10];
extern __device__ fpcomplex *Amps_DP[10];

} // namespace GooFit
