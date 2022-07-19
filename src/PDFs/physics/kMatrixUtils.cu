#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/kMatrixUtils.h>
#include <goofit/detail/compute_inverse5.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "lineshapes/Common.h"

namespace GooFit {

// For kMatrix
__device__ auto phsp_twoBody(fptype s, fptype m0, fptype m1) -> fpcomplex {
    if(1. - POW2(m0 + m1) / s > 0)
        return fpcomplex(sqrt(1. - POW2(m0 + m1) / s), 0);
    else
        return fpcomplex(0, sqrt(-(1. - POW2(m0 + m1) / s)));
}

__device__ auto phsp_fourPi(fptype s) -> fpcomplex {
    if(s > 1)
        return phsp_twoBody(s, 2 * mPiPlus, 2 * mPiPlus);
    else if(0.00051 + -0.01933 * s + 0.13851 * s * s + -0.20840 * s * s * s + -0.29744 * s * s * s * s
                + 0.13655 * s * s * s * s * s + 1.07885 * s * s * s * s * s * s
            >= 0)
        return fpcomplex(0.00051 + -0.01933 * s + 0.13851 * s * s + -0.20840 * s * s * s + -0.29744 * s * s * s * s
                             + 0.13655 * s * s * s * s * s + 1.07885 * s * s * s * s * s * s,
                         0);
    else
        return fpcomplex(0,
                         0.00051 + -0.01933 * s + 0.13851 * s * s + -0.20840 * s * s * s + -0.29744 * s * s * s * s
                             + 0.13655 * s * s * s * s * s + 1.07885 * s * s * s * s * s * s);
}

__device__ void
getCofactor(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex temp[NCHANNELS][NCHANNELS], int p, int q, int n) {
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for(int row = 0; row < n; row++) {
        for(int col = 0; col < n; col++) {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if(row != p && col != q) {
                temp[i][j++] = A[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if(j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

// Function to get adjoint of A[N][N] in adj[N][N].
__device__ void adjoint(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex adj[NCHANNELS][NCHANNELS]) {
    if(NCHANNELS == 1) {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1;
    fpcomplex temp[NCHANNELS][NCHANNELS];

    for(int i = 0; i < NCHANNELS; i++) {
        for(int j = 0; j < NCHANNELS; j++) {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, NCHANNELS);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = fptype(sign) * (determinant<NCHANNELS - 1>(temp));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
__device__ auto inverse(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex inverse[NCHANNELS][NCHANNELS]) -> bool {
    // Find determinant of A[][]
    fpcomplex det = determinant<NCHANNELS>(A);
    if(det == fpcomplex(0, 0)) {
        printf("Singular matrix, can't find its inverse\n");
        return false;
    }

    // Find adjoint
    fpcomplex adj[NCHANNELS][NCHANNELS];
    adjoint(A, adj);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for(int i = 0; i < NCHANNELS; i++)
        for(int j = 0; j < NCHANNELS; j++)
            inverse[i][j] = adj[i][j] / det;

    return true;
}

__device__ void luDecomposition(fpcomplex A[NCHANNELS][NCHANNELS],
                                fpcomplex U[NCHANNELS][NCHANNELS],
                                fpcomplex L[NCHANNELS][NCHANNELS]) {
    for(unsigned i = 0; i < NCHANNELS; i++) {
        // Upper triangular matrix
        for(unsigned k = i; k < NCHANNELS; k++) {
            fpcomplex sum(0, 0);
            for(unsigned j = 0; j < i; j++)
                sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }

        // Lower triangular.
        for(unsigned k = i; k < NCHANNELS; k++) {
            if(i == k) {
                L[i][i] = 1;
            } else {
                fpcomplex sum(0, 0);
                for(unsigned j = 0; j < i; j++)
                    sum += L[k][j] * U[j][i];
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

__device__ bool luInverse(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex inverse[NCHANNELS][NCHANNELS]) {
    fpcomplex U[NCHANNELS][NCHANNELS]    = {0};
    fpcomplex L[NCHANNELS][NCHANNELS]    = {0};
    fpcomplex Linv[NCHANNELS][NCHANNELS] = {0};
    luDecomposition(A, U, L);

    // Compute intermediate matrix Linv.
    for(int col = 0; col < NCHANNELS; col++) {
        for(int row = 0; row < NCHANNELS; row++) {
            fpcomplex sum(0, 0);
            for(int i = 0; i < NCHANNELS; i++) {
                if(i != row) {
                    sum += L[row][i] * Linv[i][col];
                }
            }
            Linv[row][col] = ((row == col ? 1. : 0.) - sum) / L[row][row];
        }
    }

    // Calculate the inverse.
    // TODO: Can this whole calculation be done in-place?
    for(int col = 0; col < NCHANNELS; col++) {
        for(int row = NCHANNELS - 1; row >= 0; row--) {
            fpcomplex sum(0, 0);
            for(int i = 0; i < NCHANNELS; i++) {
                if(i != row) {
                    sum += U[row][i] * inverse[i][col];
                }
            }
            inverse[row][col] = (Linv[row][col] - sum) / U[row][row];
        }
    }

    return true;
}

__device__ void getPropagator(const fptype kMatrix[NCHANNELS][NCHANNELS],
                              const fpcomplex phaseSpace[NCHANNELS],
                              fpcomplex F[NCHANNELS][NCHANNELS],
                              fptype adlerTerm) {
    fpcomplex tMatrix[NCHANNELS][NCHANNELS];
    tMatrix[0][0] = fpcomplex(0, 0);

    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            tMatrix[i][j] = (i == j ? 1. : 0.) - fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
        }
    }

    luInverse(tMatrix, F);
    return;
}

} // namespace GooFit
