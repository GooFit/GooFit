#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/detail/compute_inverse5.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "../lineshapes/Common.h"

namespace GooFit {

__device__ fptype twoBodyCMmom(double rMassSq, fptype d1m, fptype d2m) {
    // For A -> B + C, calculate momentum of B and C in rest frame of A.
    // PDG 38.16.

    fptype kin1 = 1 - POW2(d1m + d2m) / rMassSq;

    kin1 = kin1 >= 0 ? sqrt(kin1) : 1;

    fptype kin2 = 1 - POW2(d1m - d2m) / rMassSq;
    kin2        = kin2 >= 0 ? sqrt(kin2) : 1;

    return 0.5 * sqrt(rMassSq) * kin1 * kin2;
}

__device__ fptype dampingFactorSquare(const fptype &cmmom, const int &spin, const fptype &mRadius) {
    fptype square = mRadius * mRadius * cmmom * cmmom;
    fptype dfsq   = 1 + square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    fptype dfsqres = dfsq + 8 + 2 * square + square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}

__device__ fptype spinFactor(unsigned int spin,
                             fptype motherMass,
                             fptype daug1Mass,
                             fptype daug2Mass,
                             fptype daug3Mass,
                             fptype m12,
                             fptype m13,
                             fptype m23,
                             unsigned int cyclic_index) {
    if(0 == spin)
        return 1; // Should not cause branching since every thread evaluates the same resonance at the same time.

    /*
    // Copied from BdkDMixDalitzAmp

    fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
    fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug3Mass));
    fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));

    fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m12 : m12));
    fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m23 : m13));
    fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

    // The above, collapsed into single tests where possible.
    fptype _mA = (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass);
    fptype _mB = (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass);
    fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));

    fptype _mAC = (PAIR_23 == cyclic_index ? m13 : m23);
    fptype _mBC = (PAIR_12 == cyclic_index ? m13 : m12);
    fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    */

    // Copied from EvtDalitzReso, with assumption that pairAng convention matches pipipi0 from EvtD0mixDalitz.
    // Again, all threads should get the same branch.
    fptype _mA  = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
    fptype _mB  = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
    fptype _mC  = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
    fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12));
    fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m12 : m13));
    fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

    fptype massFactor = 1.0 / _mAB;
    fptype sFactor    = -1;
    sFactor *= ((_mBC - _mAC) + (massFactor * (motherMass * motherMass - _mC * _mC) * (_mA * _mA - _mB * _mB)));

    if(2 == spin) {
        sFactor *= sFactor;
        fptype extraterm = ((_mAB - (2 * motherMass * motherMass) - (2 * _mC * _mC))
                            + massFactor * POW2(motherMass * motherMass - _mC * _mC));
        extraterm *= ((_mAB - (2 * _mA * _mA) - (2 * _mB * _mB)) + massFactor * POW2(_mA * _mA - _mB * _mB));
        extraterm /= 3;
        sFactor -= extraterm;
    }

    return sFactor;
}

// For kMatrix
__device__ fptype phsp_twoBody(fptype s, fptype m0, fptype m1) { return sqrt(1. - POW2(m0 + m1) / s); }

__device__ fptype phsp_fourPi(fptype s) {
    if(s > 1)
        return phsp_twoBody(s, 2 * mPiPlus, 2 * mPiPlus);
    else
        return 0.00051 + -0.01933 * s + 0.13851 * s * s + -0.20840 * s * s * s + -0.29744 * s * s * s * s
               + 0.13655 * s * s * s * s * s + 1.07885 * s * s * s * s * s * s;
}

__device__ Eigen::Array<fpcomplex, NCHANNELS, NCHANNELS>
getPropagator(const Eigen::Array<fptype, NCHANNELS, NCHANNELS> &kMatrix,
              const Eigen::Matrix<fptype, 5, 1> &phaseSpace,
              fptype adlerTerm) {
    Eigen::Array<fpcomplex, NCHANNELS, NCHANNELS> tMatrix;

    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            tMatrix(i, j) = (i == j ? 1. : 0.) - fpcomplex(0, adlerTerm) * kMatrix(i, j) * phaseSpace(j);
        }
    }

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    // Here we assume that some values are 0
    return compute_inverse5<-1,
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
    return Eigen::inverse(tMatrix);
#endif
}

} // namespace GooFit
