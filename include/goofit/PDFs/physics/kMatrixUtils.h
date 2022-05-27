#include <Eigen/Core>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#define NCHANNELS 5

namespace GooFit {
__device__ auto phsp_twoBody(fptype s, fptype m0, fptype m1) -> fpcomplex;

__device__ auto phsp_fourPi(fptype s) -> fpcomplex;

__device__ void getCofactor(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex temp[NCHANNELS][NCHANNELS], int p, int q, int n);

__device__ auto determinant(fpcomplex A[NCHANNELS][NCHANNELS], int n) -> fpcomplex;

__device__ void adjoint(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex adj[NCHANNELS][NCHANNELS]);

__device__ auto inverse(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex inverse[NCHANNELS][NCHANNELS]) -> bool;

__device__ void getPropagator(const fptype kMatrix[NCHANNELS][NCHANNELS],
                              const fpcomplex phaseSpace[NCHANNELS],
                              fpcomplex F[NCHANNELS][NCHANNELS],
                              fptype adlerTerm);



} // namespace GooFit
