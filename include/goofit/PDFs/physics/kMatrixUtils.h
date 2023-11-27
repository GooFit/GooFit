#include <Eigen/Core>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#define NCHANNELS 5

namespace GooFit {
__device__ auto phsp_twoBody(fptype s, fptype m0, fptype m1) -> fpcomplex;

__device__ auto phsp_fourPi(fptype s) -> fpcomplex;

__device__ void
getCofactor(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex temp[NCHANNELS][NCHANNELS], int p, int q, int n);

__device__ void adjoint(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex adj[NCHANNELS][NCHANNELS]);

__device__ bool inverse(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex inverse[NCHANNELS][NCHANNELS]);

__device__ void getPropagator(const fptype kMatrix[NCHANNELS][NCHANNELS],
                              const fpcomplex phaseSpace[NCHANNELS],
                              fpcomplex F[NCHANNELS][NCHANNELS],
                              fptype adlerTerm);

template <int N>
inline __device__ fpcomplex determinant(fpcomplex A[NCHANNELS][NCHANNELS]) {
    fpcomplex D = 0; // Initialize result

    fpcomplex temp[NCHANNELS][NCHANNELS]; // To store cofactors

    int sign = 1; // To store sign multiplier

    // Iterate for each element of first row
    for(int f = 0; f < N; f++) {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, N);
        D += fptype(sign) * A[0][f] * determinant<N - 1>(temp);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

template <>
inline __device__ auto determinant<1>(fpcomplex A[NCHANNELS][NCHANNELS]) -> fpcomplex {
    return A[0][0];
}

} // namespace GooFit
