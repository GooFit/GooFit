#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

#define AMP_INDICES_SIZE (600)

extern __constant__ unsigned int AmpIndices[AMP_INDICES_SIZE];
extern __device__ fpcomplex *cResSF[10];
extern __device__ fpcomplex *cResSF_TD[10];
extern __device__ fpcomplex *Amps_DP[10];

} // namespace GooFit
