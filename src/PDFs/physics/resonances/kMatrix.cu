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
    // kMatrix amplitude as described in https://arxiv.org/pdf/0804.2089.pdf, compared with AmpGen implementation

    unsigned int Mpair = pc.getConstant(0);

    // parameter index
    unsigned int idx = 0;

    // Read parameters, in the same order as they are registered at the bottom of this file
    fptype sA0      = pc.getParameter(idx++);
    fptype sA       = pc.getParameter(idx++);
    fptype s0_prod  = pc.getParameter(idx++);
    fptype s0_scatt = pc.getParameter(idx++);

    fptype fscat[NCHANNELS];
    fptype pmasses[NCHANNELS];
    fptype couplings[NCHANNELS][NCHANNELS];

    fpcomplex beta[NCHANNELS];
    fpcomplex f_prod[NCHANNELS];

    for(int i = 0; i < NCHANNELS; i++) {
        fscat[i] = pc.getParameter(idx++);
        // printf("fscat(%i) = %f\n", i, fscat[i]);
    }

    // in the next two sets of parameters the index is used two times in the same line, therefore it must be incremented
    // two times afterwards
    for(int i = 0; i < NCHANNELS; i++) {
        beta[i] = fpcomplex(pc.getParameter(idx), pc.getParameter(idx + 1));
        // printf("beta[%i] = (%f,%f)\n", i, beta[i].real(), beta[i].imag());
        idx++;
        idx++;
    }

    for(int i = 0; i < NCHANNELS; i++) {
        f_prod[i] = fpcomplex(pc.getParameter(idx), pc.getParameter(idx + 1));
        // printf("f_prod[%i] = (%f,%f)\n", i, f_prod[i].real(), f_prod[i].imag());
        idx++;
        idx++;
    }

    for(int i = 0; i < NPOLES; i++) {
        for(int j = 0; j < NPOLES; j++) {
            couplings[i][j] = pc.getParameter(idx++);
            // printf("couplings(%i,%i) = %f\n", i, j, couplings[i][j]);
        }
        pmasses[i] = pc.getParameter(idx++);
        // printf("pmasses(%i) = %f\n", i, pmasses[i]);
    }

    fptype s = (PAIR_12 == Mpair ? m12 : (PAIR_13 == Mpair ? m13 : m23));
    // fptype s = 1.53351;

    // constructKMatrix

    fptype kMatrix[NCHANNELS][NCHANNELS];

    // TODO: Make sure the order (k,i,j) is correct

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            kMatrix[i][j] = 0;
            for(int k = 0; k < 5; k++)
                kMatrix[i][j] += couplings[k][i] * couplings[k][j] / (POW2(pmasses[k]) - s);
            if(i == 0 || j == 0) // Scattering term
                kMatrix[i][j] += fscat[i + j] * (1 - s0_scatt) / (s - s0_scatt);
            // printf("kMatrix(%i,%i) = %f\n",
            //       i,
            //      j,
            //     kMatrix[i][j]
            //  );
        }
    }

    fptype adlerTerm = (1. - sA0) * (s - sA * mPiPlus * mPiPlus / 2) / (s - sA0);
    // printf("adlerTerm = %f, sA0 = %f, sA = %f, sA * mPiPlus * mPiPlus / 2 = %f \n",
    //       adlerTerm,
    //       sA0,
    //       sA,
    //       sA * mPiPlus * mPiPlus / 2);

    fpcomplex phaseSpace[NCHANNELS];
    phaseSpace[0] = phsp_twoBody(s, mPiPlus, mPiPlus);
    phaseSpace[1] = phsp_twoBody(s, mKPlus, mKPlus);
    phaseSpace[2] = phsp_fourPi(s);
    phaseSpace[3] = phsp_twoBody(s, mEta, mEta);
    phaseSpace[4] = phsp_twoBody(s, mEta, mEtap);

    fpcomplex tMatrix[NCHANNELS][NCHANNELS];
    fpcomplex I[NCHANNELS][NCHANNELS];
    // printf ("tMatrix << ");
    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            tMatrix[i][j] = (i == j ? 1. : 0.) - fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
            // printf(
            //       "(%f,%f), ",
            //      tMatrix[i][j].real(),  tMatrix[i][j].imag());
        }
        // printf("\n");
    }
    fpcomplex F[NCHANNELS][NCHANNELS];
    getPropagator(kMatrix, phaseSpace, F, adlerTerm);

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            I[i][j] = fpcomplex(0, 0);
            for(int k = 0; k < 5; k++) {
                I[i][j] = I[i][j] + F[i][k] * tMatrix[k][j];
            }
            // printf(
            //       "I(%i,%i) = (%f,%f)\n",
            //        i, j, I[i][j].real(), I[i][j].imag());
            // printf(
            //        "F(%i,%i) = (%f,%f)\n",
            //        i, j, F[i][j].real(), F[i][j].imag());
        }
    }

    // TODO: calculate out
    pc.incrementIndex(1, idx, 1, 0, 1);

    fpcomplex ret(0, 0), pole(0, 0), prod(0, 0);

    for(int pterm = 0; pterm < NPOLES; pterm++) {
        fpcomplex M = 0;
        for(int i = 0; i < NCHANNELS; i++) {
            fptype coupling = couplings[pterm][i];
            M += F[0][i] * coupling;
            // printf("F[0][%i] = (%f,%f), coupling[%i][%i] = %f\n", i, F[0][i].real(), F[0][i].imag(), i, pterm,
            // couplings[pterm][i]);
        }
        // printf("M = (%f,%f), 1/(M^2-s) = %f\n", M.real(), M.imag(), 1/(POW2(pmasses[pterm]) - s));
        pole = M / (POW2(pmasses[pterm]) - s);
        ret  = ret + beta[pterm] * pole;

        // printf("F[0][%i] = (%f,%f),  (1 - s0_prod) / (s - s0_prod) = %f\n", pterm, F[0][pterm].real(),
        // F[0][pterm].imag(),  (1 - s0_prod) / (s - s0_prod));
        prod = F[0][pterm] * (1 - s0_prod) / (s - s0_prod);
        ret  = ret + f_prod[pterm] * prod;
        // printf("pole = (%f,%f), prod = (%f,%f), return = (%f,%f) \n", pole.real(), pole.imag(), prod.real(),
        // prod.imag(), ret.real(), ret.imag() );
    }

    return ret;
} // kMatrixFunction

__device__ resonance_function_ptr ptr_to_kMatrix_res = kMatrixRes;

namespace Resonances {

kMatrix::kMatrix(std::string name,
                 Variable a_r,
                 Variable a_i,
                 Variable sA0,
                 Variable sA,
                 Variable s0_prod,
                 Variable s0_scatt,
                 std::vector<Variable> beta_r,
                 std::vector<Variable> beta_i,
                 std::vector<Variable> f_prod_r,
                 std::vector<Variable> f_prod_i,
                 std::vector<Variable> fscat,
                 std::vector<Variable> poles,
                 unsigned int L,
                 unsigned int Mpair)
    : ResonancePdf("kMatrix", name, a_r, a_i) {
    if(fscat.size() != NCHANNELS)
        throw GooFit::GeneralError("You must have {} channels in fscat, not {}", NCHANNELS, fscat.size());

    if(poles.size() != NPOLES * (NPOLES + 1))
        throw GooFit::GeneralError("You must have {}x{} channels in poles, not {}", NPOLES, NPOLES + 1, poles.size());

    registerConstant(Mpair);

    registerParameter(sA0);
    registerParameter(sA);
    registerParameter(s0_prod);
    registerParameter(s0_scatt);

    for(int i = 0; i < NCHANNELS; i++) {
        registerParameter(fscat.at(i));
    }

    for(int i = 0; i < NCHANNELS; i++) {
        registerParameter(beta_r.at(i));
        registerParameter(beta_i.at(i));
    }

    for(int i = 0; i < NCHANNELS; i++) {
        registerParameter(f_prod_r.at(i));
        registerParameter(f_prod_i.at(i));
    }

    for(int i = 0; i < NPOLES * (NPOLES + 1); i++) {
        registerParameter(poles.at(i));
    }

    registerFunction("ptr_to_kMatrix_res", ptr_to_kMatrix_res);

    initialize();
}

} // namespace Resonances
} // namespace GooFit
