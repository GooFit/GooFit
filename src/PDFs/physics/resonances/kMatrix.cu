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

__device__ auto kMatrixRes(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
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

    for(double &i : fscat) {
        i = pc.getParameter(idx++);
    }

    // in the next two sets of parameters the index is used two times in the same line, therefore it must be incremented
    // two times afterwards
    for(auto &i : beta) {
        i = fpcomplex(pc.getParameter(idx), pc.getParameter(idx + 1));
        idx++;
        idx++;
    }

    for(auto &i : f_prod) {
        i = fpcomplex(pc.getParameter(idx), pc.getParameter(idx + 1));
        idx++;
        idx++;
    }

    for(int i = 0; i < NPOLES; i++) {
        for(int j = 0; j < NPOLES; j++) {
            couplings[i][j] = pc.getParameter(idx++);
        }
        pmasses[i] = pc.getParameter(idx++);
    }

    fptype s = (PAIR_12 == Mpair ? m12 : (PAIR_13 == Mpair ? m13 : m23));

    // constructKMatrix
    fptype kMatrix[NCHANNELS][NCHANNELS];

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            kMatrix[i][j] = 0;
            for(int k = 0; k < 5; k++)
                kMatrix[i][j] += couplings[k][i] * couplings[k][j] / (POW2(pmasses[k]) - s);
            if(i == 0 || j == 0) // Scattering term
                kMatrix[i][j] += fscat[i + j] * (1 - s0_scatt) / (s - s0_scatt);
        }
    }

    fptype adlerTerm = (1. - sA0) * (s - sA * mPiPlus * mPiPlus / 2) / (s - sA0);

    fpcomplex phaseSpace[NCHANNELS];
    phaseSpace[0] = phsp_twoBody(s, mPiPlus, mPiPlus);
    phaseSpace[1] = phsp_twoBody(s, mKPlus, mKPlus);
    phaseSpace[2] = phsp_fourPi(s);
    phaseSpace[3] = phsp_twoBody(s, mEta, mEta);
    phaseSpace[4] = phsp_twoBody(s, mEta, mEtap);

    fpcomplex F[NCHANNELS][NCHANNELS];
    getPropagator(kMatrix, phaseSpace, F, adlerTerm);

    // calculates output
    pc.incrementIndex(1, idx, 1, 0, 1);

    fpcomplex ret(0, 0), pole(0, 0), prod(0, 0);

    for(int pterm = 0; pterm < NPOLES; pterm++) {
        fpcomplex M = 0;
        for(int i = 0; i < NCHANNELS; i++) {
            fptype coupling = couplings[pterm][i];
            M += F[0][i] * coupling;
        }
        pole = M / (POW2(pmasses[pterm]) - s);
        ret  = ret + beta[pterm] * pole;

        prod = F[0][pterm] * (1 - s0_prod) / (s - s0_prod);
        ret  = ret + f_prod[pterm] * prod;
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
