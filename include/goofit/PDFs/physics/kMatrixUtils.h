#include <Eigen/Core>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#define NCHANNELS 5

namespace GooFit {
__device__ auto phsp_twoBody(fptype s, fptype m0, fptype m1) -> fpcomplex;

__device__ auto phsp_fourPi(fptype s) -> fpcomplex;

__device__ void getCofactor(fptype A[NCHANNELS][NCHANNELS], fptype temp[NCHANNELS][NCHANNELS], int p, int q, int n);

__device__ auto determinant(fptype A[NCHANNELS][NCHANNELS], int n) -> fptype;

__device__ void adjoint(fptype A[NCHANNELS][NCHANNELS], fptype adj[NCHANNELS][NCHANNELS]);

__device__ auto inverse(fptype A[NCHANNELS][NCHANNELS], fptype inverse[NCHANNELS][NCHANNELS]) -> bool;

__device__ void getPropagator(const fptype kMatrix[NCHANNELS][NCHANNELS],
                              const fpcomplex phaseSpace[NCHANNELS],
                              fpcomplex F[NCHANNELS][NCHANNELS],
                              fptype adlerTerm);

} // namespace GooFit
