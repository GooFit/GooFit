#include <goofit/PDFs/physics/resonances/kMatrix.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/kMatrixUtils.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <assert.h>

//#include <Eigen/Core>
//#include <Eigen/LU>

#include "../lineshapes/Common.h"

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <goofit/detail/compute_inverse5.h>
#endif

namespace GooFit {

__device__ fpcomplex kMatrixRes(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    unsigned int pterm = pc.getConstant(0);
    bool is_pole       = pc.getConstant(1) == 1;

    printf("pterm = %i, is_pole = %d, getNumConstants() = %i \n",pterm,is_pole, pc.getNumConstants() );

    unsigned int idx = 0;
    fptype sA0       = pc.getParameter(idx++);
    fptype sA        = pc.getParameter(idx++);
    fptype s0_prod   = pc.getParameter(idx++);
    fptype s0_scatt  = pc.getParameter(idx++);

    fptype fscat[NCHANNELS];
    fptype pmasses[NPOLES];
    fptype couplings[NPOLES][NPOLES];

    for(int i = 0; i < NCHANNELS; i++) {
        fscat[i] = pc.getParameter(idx++);
    }

    for(int i = 0; i < NPOLES; i++) {
        for(int j = 0; j < NPOLES; j++)
            couplings[i][j] = pc.getParameter(idx++);
        pmasses[i] = pc.getParameter(idx++);
    }

    pc.incrementIndex(1, idx, 2, 0, 1);
    printf ("m12 = %f", m12);
    fptype s = POW2(m12);

    // constructKMatrix

    fptype kMatrix[NCHANNELS][NCHANNELS] = {0};

    // TODO: Make sure the order (k,i,j) is correct

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            for(int k = 0; k < 5; k++){
                kMatrix[i][ j] += couplings[k][ i] * couplings[k][ j] / (POW2(pmasses[k]) - s);
                printf("kMatrix(%d,%d) =  %f , pmasses(k) = %f\n", i, j, kMatrix[i][ j], pmasses[k]);
            }
            if(i == 0 || j == 0) // Scattering term
                kMatrix[i][j] += fscat[i + j] * (1 - s0_scatt) / (s - s0_scatt);
        }
    }

    fptype adlerTerm = (1. - sA0) * (s - sA * mPiPlus * mPiPlus / 2) / (s - sA0);

    fpcomplex phaseSpace[NCHANNELS];
    phaseSpace[0] = phsp_twoBody(s, mPiPlus, mPiPlus);
    phaseSpace[1] = phsp_twoBody(s, mKPlus, mKPlus);
    phaseSpace[2] = phsp_fourPi(s);
    phaseSpace[4] = phsp_twoBody(s, mEta, mEta); 
    phaseSpace[5] = phsp_twoBody(s, mEta, mEtap);
    printf("input is:              %f %f  %f\n", s, mEta, mEtap );
    fpcomplex teste = phsp_twoBody(s, mEta, mEtap);
    printf("phase space matrix is: %f %f \n", phaseSpace[0].real(), phaseSpace[0].imag() );
    printf("                       %f %f \n", phaseSpace[1].real(), phaseSpace[1].imag() );
    printf("                       %f %f \n", phaseSpace[2].real(), phaseSpace[2].imag() );
    printf("                       %f %f \n", phaseSpace[3].real(), phaseSpace[3].imag() );
    printf("                       %f %f \n", phaseSpace[4].real(), phaseSpace[4].imag() );
    printf("phsp_twoBody(s, mEta, mEtap) %f %f \n", teste.real(), teste.imag() );
    printf("adlerTerm =  %f \n", adlerTerm );

    fpcomplex F[NCHANNELS][NCHANNELS],tMatrix[NCHANNELS][NCHANNELS];
    //return fpcomplex(1,0);
//    getPropagator(kMatrix, phaseSpace, adlerTerm, F);
    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
           if (i==j) tMatrix[i][j] = 1. - fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
           else tMatrix[i][j] = -fpcomplex(0, adlerTerm) * kMatrix[i][j] * phaseSpace[j];
            printf("tMatrix(%i,%i) = (%f, %f) \n", i, j, tMatrix[i][j].real(), tMatrix[i][j].imag());
        }
    }
    //return phaseSpace[pterm];
    //return 1. - fpcomplex(0, adlerTerm) * kMatrix[0][pterm] * phaseSpace[pterm];

    printf("pmasses(pterm) = %f, s = %f, pterm = %d \n",pmasses[pterm], s, pterm);
    //return phaseSpace[pterm];
    //return pterm;
    assert(1==pterm);
    return phaseSpace[1];
    return kMatrix[0][pterm];
//#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    // Here we assume that some values are 0
    F[0][0] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[0][1] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[0][2] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[0][3] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[0][4] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[1][0] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[1][1] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[1][2] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[1][3] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[1][4] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[2][0] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[2][1] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[2][2] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[2][3] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[2][4] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[3][0] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[3][1] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[3][2] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[3][3] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[3][4] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[4][0] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[4][1] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[4][2] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[4][3] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);

    F[4][4] = (tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
               - tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
               + tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
               - tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
               + tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
               + tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
               - tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1])
              / (tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 + tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][0] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 - tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][0] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][3]
                 + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 - tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][1] * tMatrix[1][3] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 + tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][1] * tMatrix[1][4] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][4] * tMatrix[4][2]
                 + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][0] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][4] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][4] * tMatrix[4][2]
                 - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][4] + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][4] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][1] * tMatrix[2][4] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 + tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] - tMatrix[0][3] * tMatrix[1][4] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][1] * tMatrix[3][3] * tMatrix[4][2]
                 - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][1] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][0] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][1]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][3] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][0] * tMatrix[3][3] * tMatrix[4][2]
                 + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][3] - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][2] * tMatrix[3][3] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][0] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][1] * tMatrix[2][3] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][1] * tMatrix[4][2] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][0] * tMatrix[3][2] * tMatrix[4][1]
                 + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][0] * tMatrix[4][2] - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][1] * tMatrix[3][2] * tMatrix[4][0]
                 - tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][0] * tMatrix[4][1] + tMatrix[0][4] * tMatrix[1][3] * tMatrix[2][2] * tMatrix[3][1] * tMatrix[4][0]);
//#else
//    F = Eigen::inverse(tMatrix);
//#endif

    // TODO: calculate out

    printf("pmasses(pterm) = %f, s = %f \n",pmasses[pterm], s);
    //return fpcomplex(1,0);
    printf(is_pole ? "is_pole = true \n" : "is_pole = false \n");
    bool dummy = false;
    printf(dummy ? "A - dummy = true \n" : "A - dummy = false \n");
    if(is_pole) dummy = true;
    printf(dummy ? "B - dummy = true \n" : "B - dummy = false \n");
    for(unsigned int i = 0; i < NCHANNELS; ++i) {
        for(unsigned int j = 0; j < NCHANNELS; ++j) {
            printf("F(%i,%i) = (%f, %f) \n", i, j, F[i][j].real(), F[i][j].imag());
        }
    }
    if(is_pole) { // pole
        fpcomplex M = 0;
        printf("pmasses(pterm) = %f, s = %f \n",pmasses[pterm], s);
        for(int i = 0; i < NCHANNELS; i++) {
            fptype pole = couplings[i][pterm];
            M += F[0][i] * pole;
        }
        fpcomplex ret = M / (POW2(pmasses[pterm]) - s);
        printf("returning: (%f,%f) \n", ret.real(), ret.imag());
        printf("M = (%f,%f), pmasses(pterm) = %f, s = %f", M.real(),M.imag(),pmasses[pterm], s);
        printf("returning: (%f,%f) \n", ret.real(), ret.imag());
        //return M;
        return M / (POW2(pmasses[pterm]) - s);

    } else { // prod
        //return fpcomplex(1,0);
        //printf("s0_prod = %f, s = %f\n", s0_prod, s);
        //        printf("F(0, pterm) = (%f,%f), s0_prod = %f, s = %f\n", tMatrix(0, pterm).real(), F(0, pterm).imag(), s0_prod, s);
        printf("returning:\n");
        //return fpcomplex(1,0);
        return F[0][pterm];
        //return F[0][pterm] * (1 - s0_prod) / (s - s0_prod);
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

            printf("pterm = %i, ispole = %i \n", pterm, is_pole);

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
