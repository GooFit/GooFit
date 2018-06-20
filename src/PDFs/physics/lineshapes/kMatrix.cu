#include <goofit/PDFs/physics/lineshapes/kMatrix.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "Common.h"

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <goofit/detail/compute_inverse5.h>
#endif

#define NPOLES 5
#define NCHANNELS 5

namespace GooFit {

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

__device__ fpcomplex kMatrixFunction(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) {
    // const fptype mass  = GOOFIT_GET_PARAM(2);
    // const fptype width = GOOFIT_GET_PARAM(3);
    // const unsigned int L = GOOFIT_GET_INT(4);
    // const fptype radius = GOOFIT_GET_CONST(7);

    // const fptype pTerm = GOOFIT_GET_INT();

    unsigned int pterm = pc.getConstant(0);
    bool is_pole       = pc.getConstant(1) == 1;

    unsigned int idx = 0;
    fptype sA0       = pc.getParameter(idx++);
    fptype sA        = pc.getParameter(idx++);
    fptype s0_prod   = pc.getParameter(idx++);
    fptype s0_scatt  = pc.getParameter(idx++);

    Eigen::Array<fptype, NCHANNELS, 1> fscat;
    Eigen::Array<fptype, NPOLES, 1> pmasses;
    Eigen::Array<fptype, NPOLES, NPOLES> couplings;

    for(int i = 0; i < NCHANNELS; i++) {
        fscat(i) = pc.getParameter(idx++);
    }

    for(int i = 0; i < NPOLES; i++) {
        for(int j = 0; j < NPOLES; j++)
            couplings(i, j) = pc.getParameter(idx++);
        pmasses(i) = pc.getParameter(idx++);
    }

    fptype s = POW2(Mpair);

    // constructKMatrix

    Eigen::Array<fptype, NCHANNELS, NCHANNELS> kMatrix;
    kMatrix.setZero();

    // TODO: Make sure the order (k,i,j) is correct

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            for(int k = 0; k < 5; k++)
                kMatrix(i, j) += couplings(k, i) * couplings(k, j) / (pmasses(k) - s);
            if(i == 0 || j == 0) // Scattering term
                kMatrix(i, j) += fscat(i + j) * (1 - s0_scatt) / (s - s0_scatt);
        }
    }

    fptype adlerTerm = (1. - sA0) * (s - sA * mPiPlus * mPiPlus / 2) / (s - sA0);

    Eigen::Matrix<fptype, 5, 1> phaseSpace;
    phaseSpace << phsp_twoBody(s, mPiPlus, mPiPlus), phsp_twoBody(s, mKPlus, mKPlus), phsp_fourPi(s),
        phsp_twoBody(s, mEta, mEta), phsp_twoBody(s, mEta, mEtap);

    Eigen::Array<fpcomplex, NCHANNELS, NCHANNELS> F = getPropagator(kMatrix, phaseSpace, adlerTerm);

    // TODO: calculate out
    pc.incrementIndex(1, idx, 2, 0, 1);

    if(is_pole) { // pole
        fpcomplex M = 0;
        for(int i = 0; i < NCHANNELS; i++) {
            fptype pole = couplings(i, pterm);
            M += F(0, i) * pole;
        }
        return M / (POW2(pmasses(pterm)) - s);
    } else { // prod
        return F(0, pterm) * (1 - s0_prod) / (s - s0_prod);
    }
}

__device__ resonance_function_ptr ptr_to_kMatrix = kMatrixFunction;

Lineshapes::kMatrix::kMatrix(std::string name,
                             unsigned int pterm,
                             bool is_pole,
                             Variable sA0,
                             Variable sA,
                             Variable s0_prod,
                             Variable s0_scatt,
                             std::array<Variable, NCHANNELS> fscat,
                             std::array<Variable, NPOLES *(NPOLES + 1)> poles,
                             Variable mass,
                             Variable width,
                             unsigned int L,
                             unsigned int Mpair,
                             FF FormFac,
                             fptype radius)
    : Lineshape("kMatrix", name, L, Mpair, FormFac, radius) {
    registerConstant(pterm);
    registerConstant(is_pole ? 1 : 0);

    registerParameter(sA0);
    registerParameter(sA);
    registerParameter(s0_prod);
    registerParameter(s0_scatt);

    for(int i = 0; i < NCHANNELS; i++) {
        registerParameter(fscat.at(i));
    }

    for(int i = 0; i < NPOLES * (NPOLES + 1); i++) {
        registerParameter(poles.at(i));
    }

    registerFunction("ptr_to_kMatrix", ptr_to_kMatrix);

    initialize();
}

} // namespace GooFit
