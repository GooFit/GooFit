#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/detail/Globals.h>

namespace GooFit {

// Device-side, translation-unit constrained.  These were constant, removing const.
// The reason is that this will make it much more difficult to fetch memory, since
// it has less memory to work with limiting the constant amount.
__device__ fptype d_parameters[GOOFIT_MAXPAR];
__device__ fptype d_constants[GOOFIT_MAXPAR];
__device__ fptype d_observables[GOOFIT_MAXPAR];
__device__ fptype d_normalizations[GOOFIT_MAXPAR];

// Function-pointer related.
__device__ void *device_function_table[GOOFIT_MAXFUNC];
// Not clear why this cannot be __constant__, but it causes crashes to declare it so.

__constant__ unsigned int c_totalEvents;
__constant__ fptype c_motherMass;
__constant__ fptype c_daug1Mass;
__constant__ fptype c_daug2Mass;
__constant__ fptype c_daug3Mass;
__constant__ fptype c_meson_radius;

__device__ int dev_powi(int base, int exp) {
    int ret = 1;

    for(int i = 0; i < exp; ++i)
        ret *= base;

    return ret;
}

} // namespace GooFit
