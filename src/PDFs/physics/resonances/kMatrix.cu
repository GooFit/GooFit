#include <goofit/PDFs/physics/resonances/kMatrix.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/kMatrixUtils.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "../lineshapes/Common.h"

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <goofit/detail/compute_inverse5.h>
#endif

namespace GooFit {

__device__ fpcomplex kMatrixRes(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
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

    fptype s = 1.53351;

    // constructKMatrix

    Eigen::Array<fptype, NCHANNELS, NCHANNELS> kMatrix;
    kMatrix.setZero();

    // TODO: Make sure the order (k,i,j) is correct

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            for(int k = 0; k < 5; k++)
                kMatrix(i, j) += couplings(k, i) * couplings(k, j) / (POW2(pmasses(k)) - s);
            if(i == 0 || j == 0) // Scattering term
                kMatrix(i, j) += fscat(i + j) * (1 - s0_scatt) / (s - s0_scatt);
            for(int k = 0; k < 5; k++)
                printf("couplings(%i,%i) = %f, couplings = %f, pmasses= %f, fscat = %f, s0_scat = %f \n",
                       i,
                       j,
                       couplings(k, i),
                       couplings(i, k),
                       pmasses(k),
                       fscat(i + j),
                       s0_scatt);
        }
    }

    fptype adlerTerm = (1. - sA0) * (s - sA * mPiPlus * mPiPlus / 2) / (s - sA0);
    printf("adlerTerm = %f, sA0 = %f, sA = %f, sA * mPiPlus * mPiPlus / 2 = %f \n",
           adlerTerm,
           sA0,
           sA,
           sA * mPiPlus * mPiPlus / 2);

    Eigen::Matrix<fpcomplex, 5, 1> phaseSpace;
    phaseSpace(0, 1) = phsp_twoBody(s, mPiPlus, mPiPlus);
    phaseSpace(1, 1) = phsp_twoBody(s, mKPlus, mKPlus);
    phaseSpace(2, 1) = phsp_fourPi(s);
    phaseSpace(3, 1) = phsp_twoBody(s, mEta, mEta);
    phaseSpace(4, 1) = phsp_twoBody(s, mEta, mEtap);

    Eigen::Matrix<fpcomplex, NCHANNELS, NCHANNELS> tMatrix, I;
    printf ("tMatrix << ");
    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            tMatrix(i, j) = (i == j ? 1. : 0.) - fpcomplex(0, adlerTerm) * kMatrix(i, j) * phaseSpace(j,1);
                printf(
                         "(%f,%f), ", 
                        tMatrix(i, j).real(),  tMatrix(i, j).imag());
        }
        printf("\n");
    }
    Eigen::Matrix<fpcomplex, NCHANNELS, NCHANNELS> F;
    getPropagator(kMatrix, phaseSpace, F, adlerTerm);

    I =  F*tMatrix;
    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            printf(
                    "I(%i,%i) = (%f,%f)\n", 
                    i, j, I(i, j).real(),I(i, j).imag());
            I(i,j) = fpcomplex(0,0);
            for(int k = 0; k < 5; k++) {
                I(i,j) = I(i,j) + F(i,k)*tMatrix(k,j);
                printf(
                        "tMatrix(%i,%i) = (%f,%f), F(%i,%i ) = (%f,%f), I(%i,%i) = (%f,%f)\n", 
                        k, j, tMatrix(k, j).real(),  tMatrix(k, j).imag(), i, k , F(i, k).real(), F(i, k).imag(), i, j, I(i, j).real(), I(i, j).imag());
            }
        }
    }

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
} // kMatrixFunction

__device__ resonance_function_ptr ptr_to_kMatrix_res = kMatrixRes;

namespace Resonances {

    kMatrix::kMatrix(std::string name,
            unsigned int pterm,
            bool is_pole,
            Variable a_r,
            Variable a_i,
            Variable sA0,
            Variable sA,
            Variable s0_prod,
            Variable s0_scatt,
            std::vector<Variable> fscat,
            std::vector<Variable> poles,
            unsigned int L,
            unsigned int Mpair)
        : ResonancePdf("kMatrix", name, a_r, a_i) {
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

            registerFunction("ptr_to_kMatrix_res", ptr_to_kMatrix_res);

            initialize();
        }

} // namespace Resonances
} // namespace GooFit
