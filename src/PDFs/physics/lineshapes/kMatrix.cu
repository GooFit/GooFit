#include <goofit/PDFs/physics/lineshapes/kMatrix.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/kMatrixUtils.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "Common.h"

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <goofit/detail/compute_inverse5.h>
#endif

namespace GooFit {

__device__ fpcomplex kMatrixFunction(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) {
    // const fptype mass  = GOOFIT_GET_PARAM(2);
    // const fptype width = GOOFIT_GET_PARAM(3);
    // const unsigned int L = GOOFIT_GET_INT(4);
    // const fptype radius = GOOFIT_GET_CONST(7);

    // const fptype pTerm = GOOFIT_GET_INT();

    unsigned int pterm = pc.getConstant(1);
    bool is_pole       = pc.getConstant(2) == 1;

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

    Eigen::Matrix<fpcomplex, 5, 1> phaseSpace;
    phaseSpace << phsp_twoBody(s, mPiPlus, mPiPlus), phsp_twoBody(s, mKPlus, mKPlus), phsp_fourPi(s),
        phsp_twoBody(s, mEta, mEta), phsp_twoBody(s, mEta, mEtap);

    Eigen::Array<fpcomplex, NCHANNELS, NCHANNELS> F = getPropagator(kMatrix, phaseSpace, adlerTerm);

    // TODO: calculate out
    pc.incrementIndex(1, idx, 3, 0, 1);

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
                             std::vector<Variable> fscat,
                             std::vector<Variable> poles,
                             Variable mass,
                             Variable width,
                             unsigned int L,
                             unsigned int Mpair,
                             FF FormFac,
                             fptype radius)
    : Lineshape("kMatrix", name, L, Mpair, FormFac, radius) {
    if(fscat.size() != NCHANNELS)
        throw GooFit::GeneralError("You must have {} channels in fscat, not {}", NCHANNELS, fscat.size());

    if(poles.size() != NPOLES * (NPOLES + 1))
        throw GooFit::GeneralError("You must have {}x{} channels in poles, not {}", NPOLES, NPOLES + 1, poles.size());

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
