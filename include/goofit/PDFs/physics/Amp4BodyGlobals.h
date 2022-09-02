#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

extern __constant__ unsigned int AmpIndices[2][500];//assume maximum two Amp4BodyPDFs at any one time (should be changed at some point to be dynamically assigned)
extern __device__ fpcomplex *cResSF[10];
extern __device__ fpcomplex *cResSF_TD[10];
extern __device__ fpcomplex *Amps_DP[10];

//global parameter that counts number of Amp4Body PDFs that are created
extern unsigned int gCacheCount; //g represents global

} // namespace GooFit
