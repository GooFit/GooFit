#include <goofit/PDFs/physics/Amp4BodyGlobals.h>

namespace GooFit {

/*
Constant memory array to hold specific info for amplitude calculation.
First entries are the starting points in array, necessary, because number of Lineshapes(LS) or Spinfactors(SF) can vary
|start of each Amplitude| #Linshapes | #Spinfactors | LS-indices | SF-indices|
| 1 entry per Amplitude | 1 per Amp  | 1 per Amp    | #LS in Amp| #SF in Amp|
*/
__constant__ unsigned int AmpIndices[2][500]; //assume maximum two Amp4BodyPDFs at any one time (should be changed at some point to be dynamically assigned)

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!

__device__ fpcomplex *cResSF[10];
__device__ fpcomplex *cResSF_TD[10]; // TODO: combine with line above
__device__ fpcomplex *Amps_DP[10];

//global parameter that counts number of Amp4Body PDFs that are created
unsigned int gCacheCount; //g represents global
} // namespace GooFit
