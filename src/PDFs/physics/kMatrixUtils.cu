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
    printf("starting getCofactor with %d %d %d \n", p, q, n);
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
    printf("finished getCofactor with %d %d %d \n", p, q, n);
}





__device__ auto determinant2(fpcomplex A[NCHANNELS][NCHANNELS]) -> fpcomplex {
    
    fpcomplex D{1,0}; // Initialize result



    return D;
}



__device__ auto determinant3(fpcomplex A[NCHANNELS][NCHANNELS]) -> fpcomplex {
    

    fpcomplex first  = A[0][0]*A[1][1]*A[2][2] - A[2][0]*A[1][1]*A[0][2];
    fpcomplex second = A[0][1]*A[1][2]*A[2][0] - A[2][1]*A[1][2]*A[0][0];
    fpcomplex third  = A[0][2]*A[1][0]*A[2][1] - A[2][2]*A[1][0]*A[0][1];


    return first+second+third;
}


__device__ auto determinant4(fpcomplex A[NCHANNELS][NCHANNELS]) -> fpcomplex {
    
    fpcomplex D = 0; // Initialize result

    //  Base case : if matrix contains single element


    fpcomplex temp[NCHANNELS][NCHANNELS]; // To store cofactors


   
    int sign = 1; // To store sign multiplier
    // Iterate for each element of first row
    for(int f = 0; f < 4; f++) {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, 4);
        D += fptype(sign) * A[0][f] * determinant3(temp);
        printf("determinant4 %d %d %d \n", NCHANNELS, 4, f);
       // D += determinant4(temp, 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}


__device__ auto determinantDecider(fpcomplex A[NCHANNELS][NCHANNELS], int n) -> fpcomplex {
    if(n==4) return determinant4(A);
    if(n==3) return determinant3(A);
    if(n==2) return determinant2(A);
    if(n==1) return A[0][0];


}


/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
__device__ auto determinant(fpcomplex A[NCHANNELS][NCHANNELS], int n) -> fpcomplex {
    printf("starting determinant with %d \n", n);
    
    fpcomplex D = 0; // Initialize result
    
    //  Base case : if matrix contains single element
    if(n == 1)
        return A[0][0];

    fpcomplex temp[NCHANNELS][NCHANNELS]; // To store cofactors

    for(int i = 0; i < NCHANNELS; i++) {
        for(int j = 0; j < NCHANNELS; j++) {

            temp[i][j] = fpcomplex(0.,0.);
            printf("tempo %f %f \n", temp[i][j].real(), temp[i][j].imag());
        }
        
    }
   
    int sign = 1; // To store sign multiplier
    // Iterate for each element of first row
    for(int f = 0; f < n; f++) {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += fptype(sign) * A[0][f] * determinant4(temp);
        printf("determinant %d %d %d \n", NCHANNELS, n, f);
       // D += determinant4(temp, 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
__device__ void adjoint(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex adj[NCHANNELS][NCHANNELS]) {
    printf("starting adjoint\n");
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
            //sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = fptype(sign) * (determinant(temp, NCHANNELS - 1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
__device__ auto inverse(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex inverse[NCHANNELS][NCHANNELS]) -> bool {
    
    // Find determinant of A[][]
    fpcomplex det = determinant(A, NCHANNELS);
    if(det == fpcomplex(0, 0)) {
        printf("Singular matrix, can't find its inverse");
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

__device__ void getPropagator(const fptype kMatrix[NCHANNELS][NCHANNELS],
                              const fpcomplex phaseSpace[NCHANNELS],
                              fpcomplex F[NCHANNELS][NCHANNELS],
                              fptype adlerTerm) {
    return;
    fpcomplex tMatrix[NCHANNELS][NCHANNELS];

    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            //tMatrix[i][j] = (i == j ? 1. : 0.) - fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
            //tMatrix[i][j] =  - fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
            // printf("tMatrix(%i,%i) = (%f,%f), kMatrix(%i,%i) = %f, phaseSpace = (%f,%f) \n",
            //       i,
            //       j,
            //       tMatrix[i][j].real(),
            //       tMatrix[i][j].imag(),
            //       i,
            //       j,
            //       kMatrix[i][j],
            //       phaseSpace[j].real(),
            //       phaseSpace[j].imag());
        }
    } 
    return;

   

    /*#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    // Here we assume that some values are 0
        F = compute_inverse5<-1,
                                -1,
                                0,
                                -1,
                                -1,
                                -1,
                                -1,
                                0,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1,
                                -1>(tMatrix);
    #else
    */
    inverse(tMatrix, F);
    return;
    //#endif
    
}

} // namespace GooFit
