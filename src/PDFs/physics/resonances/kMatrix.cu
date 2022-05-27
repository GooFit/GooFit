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


__device__ auto inverse3(fpcomplex A[NCHANNELS][NCHANNELS], fpcomplex inverse[NCHANNELS][NCHANNELS]) -> bool {
    
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

__device__ void getPropagator3(const fptype kMatrix[NCHANNELS][NCHANNELS],
                              const fpcomplex phaseSpace[NCHANNELS],
                              fpcomplex F[NCHANNELS][NCHANNELS],
                              fptype adlerTerm) {

    fpcomplex tMatrix[NCHANNELS][NCHANNELS];

    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            //tMatrix[i][j] = (i == j ? 1. : 0.) - fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
            tMatrix[i][j] =  - fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
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
    inverse3(tMatrix, F);
    return;
    //#endif
    
}

__device__ auto kMatrixRes(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    // kMatrix amplitude as described in https://arxiv.org/pdf/0804.2089.pdf, compared with AmpGen implementation

    unsigned int Mpair = pc.getConstant(0);

    // parameter index
    int idx = 0;

    // Read parameters, in the same order as they are registered at the bottom of this file
    for(int i = 0; i < pc.getNumParameters(); i++) {
        printf("par %d %f \n", i, pc.getParameter(i));
    }
    fptype sA0      = pc.getParameter(idx++);
    fptype sA       = pc.getParameter(idx++);
    fptype s0_prod  = pc.getParameter(idx++);
    fptype s0_scatt = pc.getParameter(idx++);

    fptype fscat[NCHANNELS];
    fptype pmasses[NCHANNELS];
    fptype couplings[NCHANNELS][NCHANNELS];
    fptype myfscat = 0;

    fpcomplex beta[20];
    fpcomplex f_prod[NCHANNELS];
    int counter = 0;
     printf("couter %d\n",counter);
    for(int i = 0; i < NCHANNELS; i++) {
        
        fscat[i] = pc.getParameter(idx);
        printf("couter %d %f\n",counter, fscat[i]);
        idx++;
        counter++;
    }

    // in the next two sets of parameters the index is used two times in the same line, therefore it must be incremented
    // two times afterwards
    counter = 0;
     printf("couter %d\n",counter);
    for(int i = 0; i < NCHANNELS; i++) {
        beta[i] = fpcomplex(pc.getParameter(idx), pc.getParameter(idx + 1));
         printf("beta %d %f %f\n",i, beta[i].real(), beta[i].imag());
        idx++;
        idx++;
    }
    for(int i = 0; i < NCHANNELS; i++) {
        f_prod[i] = fpcomplex(pc.getParameter(idx), pc.getParameter(idx + 1));
        printf("f_prod %d %f %f\n",i, f_prod[i].real(), f_prod[i].imag());
        idx++;
        idx++;
    }


    for(int i = 0; i < NPOLES; i++) {
        for(int j = 0; j < NPOLES; j++) {
            couplings[i][j] = pc.getParameter(idx++);
            printf("couplings %d %d %f \n", i, j, couplings[i][j]);
        }
        pmasses[i] = pc.getParameter(idx++);
        printf("masses %d %f \n", i, pmasses[i]);
    }

    fptype s = (PAIR_12 == Mpair ? m12 : (PAIR_13 == Mpair ? m13 : m23));

        // constructKMatrix
    fptype kMatrix[NCHANNELS][NCHANNELS];

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            kMatrix[i][j] = 0;
            for(int k = 0; k < 5; k++)
                kMatrix[i][j] += couplings[k][i] * couplings[k][j] / (POW2(pmasses[k]) - s);
            int blub = i+j;
            //if(i == 0 || j == 0) printf("i j %d %d %d \n", i,j,blub); // Scattering term
            //    kMatrix[i][j] += fscat[i + j] * (1 - s0_scatt) / (s - s0_scatt);
            printf("kamtrix %d %d %f \n", i, j, kMatrix[i][j]);
        }
    }

    for(int i = 0; i < 5; i++) {
        kMatrix[i][0] += fscat[i] * (1 - s0_scatt) / (s - s0_scatt);
        kMatrix[0][i] += fscat[i] * (1 - s0_scatt) / (s - s0_scatt);
    }

for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            printf("kamtrix2 %d %d %f \n", i, j, kMatrix[i][j]);
        }}


    fptype adlerTerm = (1. - sA0) * (s - sA * mPiPlus * mPiPlus / 2) / (s - sA0);

    fpcomplex phaseSpace[NCHANNELS];
    phaseSpace[0] = phsp_twoBody(s, mPiPlus, mPiPlus);
    phaseSpace[1] = phsp_twoBody(s, mKPlus, mKPlus);
    phaseSpace[2] = phsp_fourPi(s);
    phaseSpace[3] = phsp_twoBody(s, mEta, mEta);
    phaseSpace[4] = phsp_twoBody(s, mEta, mEtap);
    for(int i = 0; i < NCHANNELS; i++){
        printf("phasepace %d %f %f \n", i, phaseSpace[i].real(), phaseSpace[i].imag());
    }
    
    fpcomplex F[NCHANNELS][NCHANNELS];
    getPropagator3(kMatrix, phaseSpace, F, adlerTerm);
    
    // calculates output


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



    //pc.incrementIndex();
//fpcomplex ret(kMatrix[0][0] + phaseSpace[0].real(), kMatrix[1][3]);
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
                 std::vector<Variable> & fscat,
                 std::vector<Variable> & beta_r,
                 std::vector<Variable> & beta_i,
                 std::vector<Variable> & f_prod_r,
                 std::vector<Variable> & f_prod_i,
                 std::vector<Variable> & poles,
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

    //initialize();
}

} // namespace Resonances
} // namespace GooFit
