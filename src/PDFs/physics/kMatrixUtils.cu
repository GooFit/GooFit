#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/kMatrixUtils.h>
#include <goofit/detail/compute_inverse5.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "lineshapes/Common.h"

namespace GooFit {

// For kMatrix
__device__ fpcomplex phsp_twoBody(fptype s, fptype m0, fptype m1) {
    if(1. - POW2(m0 + m1) / s > 0)
        return fpcomplex(sqrt(1. - POW2(m0 + m1) / s), 0);
    else
        return fpcomplex(0, sqrt(-(1. - POW2(m0 + m1) / s)));
}

__device__ fpcomplex phsp_fourPi(fptype s) {
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
getPropagator(const Eigen::Array<fptype, NCHANNELS, NCHANNELS> &kMatrix,
              const Eigen::Matrix<fpcomplex, 5, 1> &phaseSpace,
              Eigen::Matrix<fpcomplex, NCHANNELS, NCHANNELS> &F,
              fptype adlerTerm) {
    Eigen::Matrix<fpcomplex, NCHANNELS, NCHANNELS> tMatrix;

    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            tMatrix(i, j) = (i == j ? 1. : 0.) - fpcomplex(0, adlerTerm) * kMatrix(i, j) * phaseSpace(j);
            printf("tMatrix(%i,%i) = (%f,%f), phaseSpace = (%f,%f) \n",
                   i,
                   j,
                   tMatrix(i, j).real(),
                   tMatrix(i, j).imag(),
                   phaseSpace(j, 1).real(),
                   phaseSpace(j, 1).imag());
        }
    }

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
    F = tMatrix.inverse();
//#endif
}

} // namespace GooFit
