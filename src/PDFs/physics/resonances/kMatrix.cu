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
    unsigned int cyclic_index = pc.getConstant(1);
    fptype beta1_ampr         = 8.5 * cos(68.5 * fptype(3.14159265) / 180.);
    fptype beta1_ampi         = 8.5 * sin(68.5 * fptype(3.14159265) / 180.);
    fptype beta2_ampr         = 12.2 * cos(24.0 * fptype(3.14159265) / 180.);
    fptype beta2_ampi         = 12.2 * sin(24.0 * fptype(3.14159265) / 180.);
    fptype beta3_ampr         = 29.2 * cos(-0.1 * fptype(3.14159265) / 180.);
    fptype beta3_ampi         = 29.2 * sin(-0.1 * fptype(3.14159265) / 180.);
    fptype beta4_ampr         = 10.8 * cos(-51.9 * fptype(3.14159265) / 180.);
    fptype beta4_ampi         = 10.8 * sin(-51.9 * fptype(3.14159265) / 180.);
    fptype beta5_ampr         = 10.8 * cos(-51.9 * fptype(3.14159265) / 180.);
    fptype beta5_ampi         = 10.8 * sin(-51.9 * fptype(3.14159265) / 180.);
    fptype fprod11r           = 8.0 * cos(-126.0 * fptype(3.14159265) / 180.);
    fptype fprod11i           = 8.0 * sin(-126.0 * fptype(3.14159265) / 180.);
    fptype fprod12r           = 26.3 * cos(-152.3 * fptype(3.14159265) / 180.);
    fptype fprod12i           = 26.3 * sin(-152.3 * fptype(3.14159265) / 180.);
    fptype fprod13r           = 33.0 * cos(-92.3 * fptype(3.14159265) / 180.);
    fptype fprod13i           = 33.0 * sin(-92.3 * fptype(3.14159265) / 180.);
    fptype fprod14r           = 26.2 * cos(-121.4 * fptype(3.14159265) / 180.);
    fptype fprod14i           = 26.2 * sin(-121.4 * fptype(3.14159265) / 180.);
    fptype fprod15r           = 26.2 * cos(-121.4 * fptype(3.14159265) / 180.);
    fptype fprod15i           = 26.2 * sin(-121.4 * fptype(3.14159265) / 180.);
    fptype Spr0               = -0.07;

    // ****************************************************
    // **** Common tools for all components (beta1-5, fprod)
    // ****************************************************

    // Particle Masses (notation: h = \eta)
    /*
    fptype pionMass   = 0.13957018; // PDG: (139.57018 \pm 0.00035) MeV
    fptype pionMassSq = pionMass*pionMass;
    fptype kMass      = 0.497614;   // PDG: (493.677 \pm 0.016) MeV
    fptype hMass      = 0.547853;   // PDG: (547.853 \pm 0.024) MeV
    fptype hprimeMass = 0.95778;    // PDG: (957.78 \pm 0.06) MeV
    */
    // EvtGen
    fptype pionMass   = 0.13957;
    fptype pionMassSq = pionMass * pionMass;
    fptype kMass      = 0.493677; // using charged K value
    fptype hMass      = 0.54775;  // using PDG value
    fptype hprimeMass = 0.95778;

    // Invariant mass squared of resonant particles (pi+ pi- = m23)
    fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

    // Pole Masses:           pipi    KK    4pi    hh   hhprime
    fptype poleMassesSq[5] = {0.651, 1.2036, 1.55817, 1.21, 1.82206}; // GeV

    for(int m = 0; m <= 4; m++) {
        poleMassesSq[m] *= poleMassesSq[m];
    } // GeV^2

    // Other Fixed Single Variables
    fptype Ssc_0 = -3.92637; // GeV^2
    fptype sA_0  = -0.15;    // GeV^2
    fptype sA    = 1.0;

    // Defining symmetric f^scattering Matrix of elements
    fptype Fsc_matrix[25];
    // non-zero terms;
    // Columnwise pipi, KK, etc
    Fsc_matrix[0]  = 0.23399;
    Fsc_matrix[1]  = 0.15044;
    Fsc_matrix[5]  = 0.15044;
    Fsc_matrix[2]  = -0.20545;
    Fsc_matrix[10] = -0.20545;
    Fsc_matrix[3]  = 0.32825;
    Fsc_matrix[15] = 0.32825;
    Fsc_matrix[4]  = 0.35412;
    Fsc_matrix[20] = 0.35412;

    // all other terms are zero
    for(int i = 1; i <= 4; i++) {
        for(int j = 1; j <= 4; j++) {
            Fsc_matrix[5 * i + j] = 0.0;
        }
    }

    // Defining g^0 Matrix of elements
    fptype g0_matrix[25];
    // alpha, pipi
    g0_matrix[0] = 0.22889;
    g0_matrix[1] = 0.94128;
    g0_matrix[2] = 0.36856;
    g0_matrix[3] = 0.3365;
    g0_matrix[4] = 0.18171;

    // alpha, KK
    g0_matrix[5] = -0.55377;
    g0_matrix[6] = 0.55095;
    g0_matrix[7] = 0.23888;
    g0_matrix[8] = 0.40907;
    g0_matrix[9] = -0.17558;

    // alpha, 4pi
    g0_matrix[10] = 0.0;
    g0_matrix[11] = 0.0;
    g0_matrix[12] = 0.55639;
    g0_matrix[13] = 0.85679;
    g0_matrix[14] = -0.79658;

    // alpha, 2eta
    g0_matrix[15] = -0.39899;
    g0_matrix[16] = 0.39065;
    g0_matrix[17] = 0.1834;
    g0_matrix[18] = 0.19906;
    g0_matrix[19] = -0.00355;

    // alpha, eta etaprime
    g0_matrix[20] = -0.34639;
    g0_matrix[21] = 0.31503;
    g0_matrix[22] = 0.18681;
    g0_matrix[23] = -0.00984;
    g0_matrix[24] = 0.22358;

    fptype smallTerm = 1.;
    fptype prec      = 0.001; // was 0.001 (MW 13 June 2016)
    for(int pole = 0; pole < 5; pole++)
        if(fabs(poleMassesSq[pole] - rMassSq) < prec)
            smallTerm = poleMassesSq[pole] - rMassSq;
    // Adler zero
    fptype AdlerZero = (1.0 - sA_0) / (rMassSq - sA_0);
    AdlerZero *= (rMassSq - 0.5 * sA * pionMassSq);
    // K-matrix*i
    fpcomplex iK_matrix[25];
    for(int i = 0; i <= 24; i++) {
        iK_matrix[i] = fpcomplex(0., 0.);
    }

    fptype sum_term = 0.0, nonRes_term = 0.0;
    int irow, jcol, alpha;
    // K-Matrix is REAL and symmetric
    for(irow = 0; irow <= 4; irow++) {
        for(jcol = irow; jcol <= 4; jcol++) {
            nonRes_term = smallTerm * Fsc_matrix[5 * irow + jcol] * (1.0 - Ssc_0)
                          / (rMassSq - Ssc_0); // Background non-resonant contribution

            sum_term = 0.0;
            for(alpha = 0; alpha <= 4; alpha++) {
                if(fabs(poleMassesSq[alpha] - rMassSq) < prec)
                    sum_term += g0_matrix[5 * irow + alpha] * g0_matrix[5 * jcol + alpha];
                else
                    sum_term += smallTerm * g0_matrix[5 * irow + alpha] * g0_matrix[5 * jcol + alpha]
                                / (poleMassesSq[alpha] - rMassSq);
            }
            // Multiplication with i -> only imaginary part
            fpcomplex kmat             = fpcomplex(0, (sum_term + nonRes_term) * AdlerZero);
            iK_matrix[5 * jcol + irow] = kmat; // upper half + diagonal elements
            iK_matrix[5 * irow + jcol] = kmat;
        }
    }
    // Calculating pseudo propagator (I - iKp)^-1
    // Calculating Phase Spaces (p = \rho)
    fpcomplex rho[5];
    for(int d = 0; d <= 4; d++) {
        rho[d] = fpcomplex(0.0, 0.0);
    }

    rho[0] = rhoF(2 * pionMass, rMassSq);
    rho[1] = rhoF(2 * kMass, rMassSq);
    rho[2] = rhoFourPiF(pionMass, rMassSq);
    rho[3] = rhoF(2 * hMass, rMassSq);
    rho[4] = rhoF(hMass + hprimeMass, rMassSq);

    // Multiplying iK by diagonal -p matrix
    for(irow = 0; irow <= 4; irow++) {
        for(jcol = 0; jcol <= 4; jcol++) {
            iK_matrix[5 * irow + jcol] *= rho[irow];
            iK_matrix[5 * irow + jcol] *= -1;
        }
    }
    // Adding Identity matrix to obtain I-iKp
    for(int d = 0; d < 5; d++) {
        iK_matrix[d * 6] = fpcomplex(smallTerm * 1.0, 0.0) + iK_matrix[d * 6];
    }
    // Explicitly Defining Matrix to be Inverted
    fpcomplex n11 = iK_matrix[0];
    fpcomplex n12 = iK_matrix[1];
    fpcomplex n13 = iK_matrix[2];
    fpcomplex n14 = iK_matrix[3];
    fpcomplex n15 = iK_matrix[4];

    fpcomplex n21 = iK_matrix[5];
    fpcomplex n22 = iK_matrix[6];
    fpcomplex n23 = iK_matrix[7];
    fpcomplex n24 = iK_matrix[8];
    fpcomplex n25 = iK_matrix[9];

    fpcomplex n31 = iK_matrix[10];
    fpcomplex n32 = iK_matrix[11];
    fpcomplex n33 = iK_matrix[12];
    fpcomplex n34 = iK_matrix[13];
    fpcomplex n35 = iK_matrix[14];

    fpcomplex n41 = iK_matrix[15];
    fpcomplex n42 = iK_matrix[16];
    fpcomplex n43 = iK_matrix[17];
    fpcomplex n44 = iK_matrix[18];
    fpcomplex n45 = iK_matrix[19];

    fpcomplex n51 = iK_matrix[20];
    fpcomplex n52 = iK_matrix[21];
    fpcomplex n53 = iK_matrix[22];
    fpcomplex n54 = iK_matrix[23];
    fpcomplex n55 = iK_matrix[24];
    // Computing elements of the first row of (I-iKp)^-1
    // Formulae for inverted matrix elements obtained from Maple

    fpcomplex inv11(0.0, 0.0);
    fpcomplex inv12(0.0, 0.0);
    fpcomplex inv13(0.0, 0.0);
    fpcomplex inv14(0.0, 0.0);
    fpcomplex inv15(0.0, 0.0);

    inv11 += (n25 * n34 * n43 * n52 - n24 * n35 * n43 * n52 - n25 * n33 * n44 * n52 + n23 * n35 * n44 * n52
              + n24 * n33 * n45 * n52 - n23 * n34 * n45 * n52 - n25 * n34 * n42 * n53 + n24 * n35 * n42 * n53
              + n25 * n32 * n44 * n53 - n22 * n35 * n44 * n53 - n24 * n32 * n45 * n53 + n22 * n34 * n45 * n53
              + n25 * n33 * n42 * n54 - n23 * n35 * n42 * n54 - n25 * n32 * n43 * n54 + n22 * n35 * n43 * n54
              + n23 * n32 * n45 * n54 - n22 * n33 * n45 * n54 - n24 * n33 * n42 * n55 + n23 * n34 * n42 * n55
              + n24 * n32 * n43 * n55 - n22 * n34 * n43 * n55 - n23 * n32 * n44 * n55 + n22 * n33 * n44 * n55)
             / (n15 * n24 * n33 * n42 * n51 - n14 * n25 * n33 * n42 * n51 - n15 * n23 * n34 * n42 * n51
                + n13 * n25 * n34 * n42 * n51 + n14 * n23 * n35 * n42 * n51 - n13 * n24 * n35 * n42 * n51
                - n15 * n24 * n32 * n43 * n51 + n14 * n25 * n32 * n43 * n51 + n15 * n22 * n34 * n43 * n51
                - n12 * n25 * n34 * n43 * n51 - n14 * n22 * n35 * n43 * n51 + n12 * n24 * n35 * n43 * n51
                + n15 * n23 * n32 * n44 * n51 - n13 * n25 * n32 * n44 * n51 - n15 * n22 * n33 * n44 * n51
                + n12 * n25 * n33 * n44 * n51 + n13 * n22 * n35 * n44 * n51 - n12 * n23 * n35 * n44 * n51
                - n14 * n23 * n32 * n45 * n51 + n13 * n24 * n32 * n45 * n51 + n14 * n22 * n33 * n45 * n51
                - n12 * n24 * n33 * n45 * n51 - n13 * n22 * n34 * n45 * n51 + n12 * n23 * n34 * n45 * n51
                - n15 * n24 * n33 * n41 * n52 + n14 * n25 * n33 * n41 * n52 + n15 * n23 * n34 * n41 * n52
                - n13 * n25 * n34 * n41 * n52 - n14 * n23 * n35 * n41 * n52 + n13 * n24 * n35 * n41 * n52
                + n15 * n24 * n31 * n43 * n52 - n14 * n25 * n31 * n43 * n52 - n15 * n21 * n34 * n43 * n52
                + n11 * n25 * n34 * n43 * n52 + n14 * n21 * n35 * n43 * n52 - n11 * n24 * n35 * n43 * n52
                - n15 * n23 * n31 * n44 * n52 + n13 * n25 * n31 * n44 * n52 + n15 * n21 * n33 * n44 * n52
                - n11 * n25 * n33 * n44 * n52 - n13 * n21 * n35 * n44 * n52 + n11 * n23 * n35 * n44 * n52
                + n14 * n23 * n31 * n45 * n52 - n13 * n24 * n31 * n45 * n52 - n14 * n21 * n33 * n45 * n52
                + n11 * n24 * n33 * n45 * n52 + n13 * n21 * n34 * n45 * n52 - n11 * n23 * n34 * n45 * n52
                + n15 * n24 * n32 * n41 * n53 - n14 * n25 * n32 * n41 * n53 - n15 * n22 * n34 * n41 * n53
                + n12 * n25 * n34 * n41 * n53 + n14 * n22 * n35 * n41 * n53 - n12 * n24 * n35 * n41 * n53
                - n15 * n24 * n31 * n42 * n53 + n14 * n25 * n31 * n42 * n53 + n15 * n21 * n34 * n42 * n53
                - n11 * n25 * n34 * n42 * n53 - n14 * n21 * n35 * n42 * n53 + n11 * n24 * n35 * n42 * n53
                + n15 * n22 * n31 * n44 * n53 - n12 * n25 * n31 * n44 * n53 - n15 * n21 * n32 * n44 * n53
                + n11 * n25 * n32 * n44 * n53 + n12 * n21 * n35 * n44 * n53 - n11 * n22 * n35 * n44 * n53
                - n14 * n22 * n31 * n45 * n53 + n12 * n24 * n31 * n45 * n53 + n14 * n21 * n32 * n45 * n53
                - n11 * n24 * n32 * n45 * n53 - n12 * n21 * n34 * n45 * n53 + n11 * n22 * n34 * n45 * n53
                - n15 * n23 * n32 * n41 * n54 + n13 * n25 * n32 * n41 * n54 + n15 * n22 * n33 * n41 * n54
                - n12 * n25 * n33 * n41 * n54 - n13 * n22 * n35 * n41 * n54 + n12 * n23 * n35 * n41 * n54
                + n15 * n23 * n31 * n42 * n54 - n13 * n25 * n31 * n42 * n54 - n15 * n21 * n33 * n42 * n54
                + n11 * n25 * n33 * n42 * n54 + n13 * n21 * n35 * n42 * n54 - n11 * n23 * n35 * n42 * n54
                - n15 * n22 * n31 * n43 * n54 + n12 * n25 * n31 * n43 * n54 + n15 * n21 * n32 * n43 * n54
                - n11 * n25 * n32 * n43 * n54 - n12 * n21 * n35 * n43 * n54 + n11 * n22 * n35 * n43 * n54
                + n13 * n22 * n31 * n45 * n54 - n12 * n23 * n31 * n45 * n54 - n13 * n21 * n32 * n45 * n54
                + n11 * n23 * n32 * n45 * n54 + n12 * n21 * n33 * n45 * n54 - n11 * n22 * n33 * n45 * n54
                + n14 * n23 * n32 * n41 * n55 - n13 * n24 * n32 * n41 * n55 - n14 * n22 * n33 * n41 * n55
                + n12 * n24 * n33 * n41 * n55 + n13 * n22 * n34 * n41 * n55 - n12 * n23 * n34 * n41 * n55
                - n14 * n23 * n31 * n42 * n55 + n13 * n24 * n31 * n42 * n55 + n14 * n21 * n33 * n42 * n55
                - n11 * n24 * n33 * n42 * n55 - n13 * n21 * n34 * n42 * n55 + n11 * n23 * n34 * n42 * n55
                + n14 * n22 * n31 * n43 * n55 - n12 * n24 * n31 * n43 * n55 - n14 * n21 * n32 * n43 * n55
                + n11 * n24 * n32 * n43 * n55 + n12 * n21 * n34 * n43 * n55 - n11 * n22 * n34 * n43 * n55
                - n13 * n22 * n31 * n44 * n55 + n12 * n23 * n31 * n44 * n55 + n13 * n21 * n32 * n44 * n55
                - n11 * n23 * n32 * n44 * n55 - n12 * n21 * n33 * n44 * n55 + n11 * n22 * n33 * n44 * n55);

    inv12 += (n24 * n35 * n43 * n51 - n25 * n34 * n43 * n51 + n25 * n33 * n44 * n51 - n23 * n35 * n44 * n51
              - n24 * n33 * n45 * n51 + n23 * n34 * n45 * n51 + n25 * n34 * n41 * n53 - n24 * n35 * n41 * n53
              - n25 * n31 * n44 * n53 + n21 * n35 * n44 * n53 + n24 * n31 * n45 * n53 - n21 * n34 * n45 * n53
              - n25 * n33 * n41 * n54 + n23 * n35 * n41 * n54 + n25 * n31 * n43 * n54 - n21 * n35 * n43 * n54
              - n23 * n31 * n45 * n54 + n21 * n33 * n45 * n54 + n24 * n33 * n41 * n55 - n23 * n34 * n41 * n55
              - n24 * n31 * n43 * n55 + n21 * n34 * n43 * n55 + n23 * n31 * n44 * n55 - n21 * n33 * n44 * n55)
             / (n15 * n24 * n33 * n42 * n51 - n14 * n25 * n33 * n42 * n51 - n15 * n23 * n34 * n42 * n51
                + n13 * n25 * n34 * n42 * n51 + n14 * n23 * n35 * n42 * n51 - n13 * n24 * n35 * n42 * n51
                - n15 * n24 * n32 * n43 * n51 + n14 * n25 * n32 * n43 * n51 + n15 * n22 * n34 * n43 * n51
                - n12 * n25 * n34 * n43 * n51 - n14 * n22 * n35 * n43 * n51 + n12 * n24 * n35 * n43 * n51
                + n15 * n23 * n32 * n44 * n51 - n13 * n25 * n32 * n44 * n51 - n15 * n22 * n33 * n44 * n51
                + n12 * n25 * n33 * n44 * n51 + n13 * n22 * n35 * n44 * n51 - n12 * n23 * n35 * n44 * n51
                - n14 * n23 * n32 * n45 * n51 + n13 * n24 * n32 * n45 * n51 + n14 * n22 * n33 * n45 * n51
                - n12 * n24 * n33 * n45 * n51 - n13 * n22 * n34 * n45 * n51 + n12 * n23 * n34 * n45 * n51
                - n15 * n24 * n33 * n41 * n52 + n14 * n25 * n33 * n41 * n52 + n15 * n23 * n34 * n41 * n52
                - n13 * n25 * n34 * n41 * n52 - n14 * n23 * n35 * n41 * n52 + n13 * n24 * n35 * n41 * n52
                + n15 * n24 * n31 * n43 * n52 - n14 * n25 * n31 * n43 * n52 - n15 * n21 * n34 * n43 * n52
                + n11 * n25 * n34 * n43 * n52 + n14 * n21 * n35 * n43 * n52 - n11 * n24 * n35 * n43 * n52
                - n15 * n23 * n31 * n44 * n52 + n13 * n25 * n31 * n44 * n52 + n15 * n21 * n33 * n44 * n52
                - n11 * n25 * n33 * n44 * n52 - n13 * n21 * n35 * n44 * n52 + n11 * n23 * n35 * n44 * n52
                + n14 * n23 * n31 * n45 * n52 - n13 * n24 * n31 * n45 * n52 - n14 * n21 * n33 * n45 * n52
                + n11 * n24 * n33 * n45 * n52 + n13 * n21 * n34 * n45 * n52 - n11 * n23 * n34 * n45 * n52
                + n15 * n24 * n32 * n41 * n53 - n14 * n25 * n32 * n41 * n53 - n15 * n22 * n34 * n41 * n53
                + n12 * n25 * n34 * n41 * n53 + n14 * n22 * n35 * n41 * n53 - n12 * n24 * n35 * n41 * n53
                - n15 * n24 * n31 * n42 * n53 + n14 * n25 * n31 * n42 * n53 + n15 * n21 * n34 * n42 * n53
                - n11 * n25 * n34 * n42 * n53 - n14 * n21 * n35 * n42 * n53 + n11 * n24 * n35 * n42 * n53
                + n15 * n22 * n31 * n44 * n53 - n12 * n25 * n31 * n44 * n53 - n15 * n21 * n32 * n44 * n53
                + n11 * n25 * n32 * n44 * n53 + n12 * n21 * n35 * n44 * n53 - n11 * n22 * n35 * n44 * n53
                - n14 * n22 * n31 * n45 * n53 + n12 * n24 * n31 * n45 * n53 + n14 * n21 * n32 * n45 * n53
                - n11 * n24 * n32 * n45 * n53 - n12 * n21 * n34 * n45 * n53 + n11 * n22 * n34 * n45 * n53
                - n15 * n23 * n32 * n41 * n54 + n13 * n25 * n32 * n41 * n54 + n15 * n22 * n33 * n41 * n54
                - n12 * n25 * n33 * n41 * n54 - n13 * n22 * n35 * n41 * n54 + n12 * n23 * n35 * n41 * n54
                + n15 * n23 * n31 * n42 * n54 - n13 * n25 * n31 * n42 * n54 - n15 * n21 * n33 * n42 * n54
                + n11 * n25 * n33 * n42 * n54 + n13 * n21 * n35 * n42 * n54 - n11 * n23 * n35 * n42 * n54
                - n15 * n22 * n31 * n43 * n54 + n12 * n25 * n31 * n43 * n54 + n15 * n21 * n32 * n43 * n54
                - n11 * n25 * n32 * n43 * n54 - n12 * n21 * n35 * n43 * n54 + n11 * n22 * n35 * n43 * n54
                + n13 * n22 * n31 * n45 * n54 - n12 * n23 * n31 * n45 * n54 - n13 * n21 * n32 * n45 * n54
                + n11 * n23 * n32 * n45 * n54 + n12 * n21 * n33 * n45 * n54 - n11 * n22 * n33 * n45 * n54
                + n14 * n23 * n32 * n41 * n55 - n13 * n24 * n32 * n41 * n55 - n14 * n22 * n33 * n41 * n55
                + n12 * n24 * n33 * n41 * n55 + n13 * n22 * n34 * n41 * n55 - n12 * n23 * n34 * n41 * n55
                - n14 * n23 * n31 * n42 * n55 + n13 * n24 * n31 * n42 * n55 + n14 * n21 * n33 * n42 * n55
                - n11 * n24 * n33 * n42 * n55 - n13 * n21 * n34 * n42 * n55 + n11 * n23 * n34 * n42 * n55
                + n14 * n22 * n31 * n43 * n55 - n12 * n24 * n31 * n43 * n55 - n14 * n21 * n32 * n43 * n55
                + n11 * n24 * n32 * n43 * n55 + n12 * n21 * n34 * n43 * n55 - n11 * n22 * n34 * n43 * n55
                - n13 * n22 * n31 * n44 * n55 + n12 * n23 * n31 * n44 * n55 + n13 * n21 * n32 * n44 * n55
                - n11 * n23 * n32 * n44 * n55 - n12 * n21 * n33 * n44 * n55 + n11 * n22 * n33 * n44 * n55);

    inv13 += (n25 * n34 * n42 * n51 - n24 * n35 * n42 * n51 - n25 * n32 * n44 * n51 + n22 * n35 * n44 * n51
              + n24 * n32 * n45 * n51 - n22 * n34 * n45 * n51 - n25 * n34 * n41 * n52 + n24 * n35 * n41 * n52
              + n25 * n31 * n44 * n52 - n21 * n35 * n44 * n52 - n24 * n31 * n45 * n52 + n21 * n34 * n45 * n52
              + n25 * n32 * n41 * n54 - n22 * n35 * n41 * n54 - n25 * n31 * n42 * n54 + n21 * n35 * n42 * n54
              + n22 * n31 * n45 * n54 - n21 * n32 * n45 * n54 - n24 * n32 * n41 * n55 + n22 * n34 * n41 * n55
              + n24 * n31 * n42 * n55 - n21 * n34 * n42 * n55 - n22 * n31 * n44 * n55 + n21 * n32 * n44 * n55)
             / (n15 * n24 * n33 * n42 * n51 - n14 * n25 * n33 * n42 * n51 - n15 * n23 * n34 * n42 * n51
                + n13 * n25 * n34 * n42 * n51 + n14 * n23 * n35 * n42 * n51 - n13 * n24 * n35 * n42 * n51
                - n15 * n24 * n32 * n43 * n51 + n14 * n25 * n32 * n43 * n51 + n15 * n22 * n34 * n43 * n51
                - n12 * n25 * n34 * n43 * n51 - n14 * n22 * n35 * n43 * n51 + n12 * n24 * n35 * n43 * n51
                + n15 * n23 * n32 * n44 * n51 - n13 * n25 * n32 * n44 * n51 - n15 * n22 * n33 * n44 * n51
                + n12 * n25 * n33 * n44 * n51 + n13 * n22 * n35 * n44 * n51 - n12 * n23 * n35 * n44 * n51
                - n14 * n23 * n32 * n45 * n51 + n13 * n24 * n32 * n45 * n51 + n14 * n22 * n33 * n45 * n51
                - n12 * n24 * n33 * n45 * n51 - n13 * n22 * n34 * n45 * n51 + n12 * n23 * n34 * n45 * n51
                - n15 * n24 * n33 * n41 * n52 + n14 * n25 * n33 * n41 * n52 + n15 * n23 * n34 * n41 * n52
                - n13 * n25 * n34 * n41 * n52 - n14 * n23 * n35 * n41 * n52 + n13 * n24 * n35 * n41 * n52
                + n15 * n24 * n31 * n43 * n52 - n14 * n25 * n31 * n43 * n52 - n15 * n21 * n34 * n43 * n52
                + n11 * n25 * n34 * n43 * n52 + n14 * n21 * n35 * n43 * n52 - n11 * n24 * n35 * n43 * n52
                - n15 * n23 * n31 * n44 * n52 + n13 * n25 * n31 * n44 * n52 + n15 * n21 * n33 * n44 * n52
                - n11 * n25 * n33 * n44 * n52 - n13 * n21 * n35 * n44 * n52 + n11 * n23 * n35 * n44 * n52
                + n14 * n23 * n31 * n45 * n52 - n13 * n24 * n31 * n45 * n52 - n14 * n21 * n33 * n45 * n52
                + n11 * n24 * n33 * n45 * n52 + n13 * n21 * n34 * n45 * n52 - n11 * n23 * n34 * n45 * n52
                + n15 * n24 * n32 * n41 * n53 - n14 * n25 * n32 * n41 * n53 - n15 * n22 * n34 * n41 * n53
                + n12 * n25 * n34 * n41 * n53 + n14 * n22 * n35 * n41 * n53 - n12 * n24 * n35 * n41 * n53
                - n15 * n24 * n31 * n42 * n53 + n14 * n25 * n31 * n42 * n53 + n15 * n21 * n34 * n42 * n53
                - n11 * n25 * n34 * n42 * n53 - n14 * n21 * n35 * n42 * n53 + n11 * n24 * n35 * n42 * n53
                + n15 * n22 * n31 * n44 * n53 - n12 * n25 * n31 * n44 * n53 - n15 * n21 * n32 * n44 * n53
                + n11 * n25 * n32 * n44 * n53 + n12 * n21 * n35 * n44 * n53 - n11 * n22 * n35 * n44 * n53
                - n14 * n22 * n31 * n45 * n53 + n12 * n24 * n31 * n45 * n53 + n14 * n21 * n32 * n45 * n53
                - n11 * n24 * n32 * n45 * n53 - n12 * n21 * n34 * n45 * n53 + n11 * n22 * n34 * n45 * n53
                - n15 * n23 * n32 * n41 * n54 + n13 * n25 * n32 * n41 * n54 + n15 * n22 * n33 * n41 * n54
                - n12 * n25 * n33 * n41 * n54 - n13 * n22 * n35 * n41 * n54 + n12 * n23 * n35 * n41 * n54
                + n15 * n23 * n31 * n42 * n54 - n13 * n25 * n31 * n42 * n54 - n15 * n21 * n33 * n42 * n54
                + n11 * n25 * n33 * n42 * n54 + n13 * n21 * n35 * n42 * n54 - n11 * n23 * n35 * n42 * n54
                - n15 * n22 * n31 * n43 * n54 + n12 * n25 * n31 * n43 * n54 + n15 * n21 * n32 * n43 * n54
                - n11 * n25 * n32 * n43 * n54 - n12 * n21 * n35 * n43 * n54 + n11 * n22 * n35 * n43 * n54
                + n13 * n22 * n31 * n45 * n54 - n12 * n23 * n31 * n45 * n54 - n13 * n21 * n32 * n45 * n54
                + n11 * n23 * n32 * n45 * n54 + n12 * n21 * n33 * n45 * n54 - n11 * n22 * n33 * n45 * n54
                + n14 * n23 * n32 * n41 * n55 - n13 * n24 * n32 * n41 * n55 - n14 * n22 * n33 * n41 * n55
                + n12 * n24 * n33 * n41 * n55 + n13 * n22 * n34 * n41 * n55 - n12 * n23 * n34 * n41 * n55
                - n14 * n23 * n31 * n42 * n55 + n13 * n24 * n31 * n42 * n55 + n14 * n21 * n33 * n42 * n55
                - n11 * n24 * n33 * n42 * n55 - n13 * n21 * n34 * n42 * n55 + n11 * n23 * n34 * n42 * n55
                + n14 * n22 * n31 * n43 * n55 - n12 * n24 * n31 * n43 * n55 - n14 * n21 * n32 * n43 * n55
                + n11 * n24 * n32 * n43 * n55 + n12 * n21 * n34 * n43 * n55 - n11 * n22 * n34 * n43 * n55
                - n13 * n22 * n31 * n44 * n55 + n12 * n23 * n31 * n44 * n55 + n13 * n21 * n32 * n44 * n55
                - n11 * n23 * n32 * n44 * n55 - n12 * n21 * n33 * n44 * n55 + n11 * n22 * n33 * n44 * n55);

    inv14 += (n23 * n35 * n42 * n51 - n25 * n33 * n42 * n51 + n25 * n32 * n43 * n51 - n22 * n35 * n43 * n51
              - n23 * n32 * n45 * n51 + n22 * n33 * n45 * n51 + n25 * n33 * n41 * n52 - n23 * n35 * n41 * n52
              - n25 * n31 * n43 * n52 + n21 * n35 * n43 * n52 + n23 * n31 * n45 * n52 - n21 * n33 * n45 * n52
              - n25 * n32 * n41 * n53 + n22 * n35 * n41 * n53 + n25 * n31 * n42 * n53 - n21 * n35 * n42 * n53
              - n22 * n31 * n45 * n53 + n21 * n32 * n45 * n53 + n23 * n32 * n41 * n55 - n22 * n33 * n41 * n55
              - n23 * n31 * n42 * n55 + n21 * n33 * n42 * n55 + n22 * n31 * n43 * n55 - n21 * n32 * n43 * n55)
             / (n15 * n24 * n33 * n42 * n51 - n14 * n25 * n33 * n42 * n51 - n15 * n23 * n34 * n42 * n51
                + n13 * n25 * n34 * n42 * n51 + n14 * n23 * n35 * n42 * n51 - n13 * n24 * n35 * n42 * n51
                - n15 * n24 * n32 * n43 * n51 + n14 * n25 * n32 * n43 * n51 + n15 * n22 * n34 * n43 * n51
                - n12 * n25 * n34 * n43 * n51 - n14 * n22 * n35 * n43 * n51 + n12 * n24 * n35 * n43 * n51
                + n15 * n23 * n32 * n44 * n51 - n13 * n25 * n32 * n44 * n51 - n15 * n22 * n33 * n44 * n51
                + n12 * n25 * n33 * n44 * n51 + n13 * n22 * n35 * n44 * n51 - n12 * n23 * n35 * n44 * n51
                - n14 * n23 * n32 * n45 * n51 + n13 * n24 * n32 * n45 * n51 + n14 * n22 * n33 * n45 * n51
                - n12 * n24 * n33 * n45 * n51 - n13 * n22 * n34 * n45 * n51 + n12 * n23 * n34 * n45 * n51
                - n15 * n24 * n33 * n41 * n52 + n14 * n25 * n33 * n41 * n52 + n15 * n23 * n34 * n41 * n52
                - n13 * n25 * n34 * n41 * n52 - n14 * n23 * n35 * n41 * n52 + n13 * n24 * n35 * n41 * n52
                + n15 * n24 * n31 * n43 * n52 - n14 * n25 * n31 * n43 * n52 - n15 * n21 * n34 * n43 * n52
                + n11 * n25 * n34 * n43 * n52 + n14 * n21 * n35 * n43 * n52 - n11 * n24 * n35 * n43 * n52
                - n15 * n23 * n31 * n44 * n52 + n13 * n25 * n31 * n44 * n52 + n15 * n21 * n33 * n44 * n52
                - n11 * n25 * n33 * n44 * n52 - n13 * n21 * n35 * n44 * n52 + n11 * n23 * n35 * n44 * n52
                + n14 * n23 * n31 * n45 * n52 - n13 * n24 * n31 * n45 * n52 - n14 * n21 * n33 * n45 * n52
                + n11 * n24 * n33 * n45 * n52 + n13 * n21 * n34 * n45 * n52 - n11 * n23 * n34 * n45 * n52
                + n15 * n24 * n32 * n41 * n53 - n14 * n25 * n32 * n41 * n53 - n15 * n22 * n34 * n41 * n53
                + n12 * n25 * n34 * n41 * n53 + n14 * n22 * n35 * n41 * n53 - n12 * n24 * n35 * n41 * n53
                - n15 * n24 * n31 * n42 * n53 + n14 * n25 * n31 * n42 * n53 + n15 * n21 * n34 * n42 * n53
                - n11 * n25 * n34 * n42 * n53 - n14 * n21 * n35 * n42 * n53 + n11 * n24 * n35 * n42 * n53
                + n15 * n22 * n31 * n44 * n53 - n12 * n25 * n31 * n44 * n53 - n15 * n21 * n32 * n44 * n53
                + n11 * n25 * n32 * n44 * n53 + n12 * n21 * n35 * n44 * n53 - n11 * n22 * n35 * n44 * n53
                - n14 * n22 * n31 * n45 * n53 + n12 * n24 * n31 * n45 * n53 + n14 * n21 * n32 * n45 * n53
                - n11 * n24 * n32 * n45 * n53 - n12 * n21 * n34 * n45 * n53 + n11 * n22 * n34 * n45 * n53
                - n15 * n23 * n32 * n41 * n54 + n13 * n25 * n32 * n41 * n54 + n15 * n22 * n33 * n41 * n54
                - n12 * n25 * n33 * n41 * n54 - n13 * n22 * n35 * n41 * n54 + n12 * n23 * n35 * n41 * n54
                + n15 * n23 * n31 * n42 * n54 - n13 * n25 * n31 * n42 * n54 - n15 * n21 * n33 * n42 * n54
                + n11 * n25 * n33 * n42 * n54 + n13 * n21 * n35 * n42 * n54 - n11 * n23 * n35 * n42 * n54
                - n15 * n22 * n31 * n43 * n54 + n12 * n25 * n31 * n43 * n54 + n15 * n21 * n32 * n43 * n54
                - n11 * n25 * n32 * n43 * n54 - n12 * n21 * n35 * n43 * n54 + n11 * n22 * n35 * n43 * n54
                + n13 * n22 * n31 * n45 * n54 - n12 * n23 * n31 * n45 * n54 - n13 * n21 * n32 * n45 * n54
                + n11 * n23 * n32 * n45 * n54 + n12 * n21 * n33 * n45 * n54 - n11 * n22 * n33 * n45 * n54
                + n14 * n23 * n32 * n41 * n55 - n13 * n24 * n32 * n41 * n55 - n14 * n22 * n33 * n41 * n55
                + n12 * n24 * n33 * n41 * n55 + n13 * n22 * n34 * n41 * n55 - n12 * n23 * n34 * n41 * n55
                - n14 * n23 * n31 * n42 * n55 + n13 * n24 * n31 * n42 * n55 + n14 * n21 * n33 * n42 * n55
                - n11 * n24 * n33 * n42 * n55 - n13 * n21 * n34 * n42 * n55 + n11 * n23 * n34 * n42 * n55
                + n14 * n22 * n31 * n43 * n55 - n12 * n24 * n31 * n43 * n55 - n14 * n21 * n32 * n43 * n55
                + n11 * n24 * n32 * n43 * n55 + n12 * n21 * n34 * n43 * n55 - n11 * n22 * n34 * n43 * n55
                - n13 * n22 * n31 * n44 * n55 + n12 * n23 * n31 * n44 * n55 + n13 * n21 * n32 * n44 * n55
                - n11 * n23 * n32 * n44 * n55 - n12 * n21 * n33 * n44 * n55 + n11 * n22 * n33 * n44 * n55);

    inv15 += (n24 * n33 * n42 * n51 - n23 * n34 * n42 * n51 - n24 * n32 * n43 * n51 + n22 * n34 * n43 * n51
              + n23 * n32 * n44 * n51 - n22 * n33 * n44 * n51 - n24 * n33 * n41 * n52 + n23 * n34 * n41 * n52
              + n24 * n31 * n43 * n52 - n21 * n34 * n43 * n52 - n23 * n31 * n44 * n52 + n21 * n33 * n44 * n52
              + n24 * n32 * n41 * n53 - n22 * n34 * n41 * n53 - n24 * n31 * n42 * n53 + n21 * n34 * n42 * n53
              + n22 * n31 * n44 * n53 - n21 * n32 * n44 * n53 - n23 * n32 * n41 * n54 + n22 * n33 * n41 * n54
              + n23 * n31 * n42 * n54 - n21 * n33 * n42 * n54 - n22 * n31 * n43 * n54 + n21 * n32 * n43 * n54)
             / (n15 * n24 * n33 * n42 * n51 - n14 * n25 * n33 * n42 * n51 - n15 * n23 * n34 * n42 * n51
                + n13 * n25 * n34 * n42 * n51 + n14 * n23 * n35 * n42 * n51 - n13 * n24 * n35 * n42 * n51
                - n15 * n24 * n32 * n43 * n51 + n14 * n25 * n32 * n43 * n51 + n15 * n22 * n34 * n43 * n51
                - n12 * n25 * n34 * n43 * n51 - n14 * n22 * n35 * n43 * n51 + n12 * n24 * n35 * n43 * n51
                + n15 * n23 * n32 * n44 * n51 - n13 * n25 * n32 * n44 * n51 - n15 * n22 * n33 * n44 * n51
                + n12 * n25 * n33 * n44 * n51 + n13 * n22 * n35 * n44 * n51 - n12 * n23 * n35 * n44 * n51
                - n14 * n23 * n32 * n45 * n51 + n13 * n24 * n32 * n45 * n51 + n14 * n22 * n33 * n45 * n51
                - n12 * n24 * n33 * n45 * n51 - n13 * n22 * n34 * n45 * n51 + n12 * n23 * n34 * n45 * n51
                - n15 * n24 * n33 * n41 * n52 + n14 * n25 * n33 * n41 * n52 + n15 * n23 * n34 * n41 * n52
                - n13 * n25 * n34 * n41 * n52 - n14 * n23 * n35 * n41 * n52 + n13 * n24 * n35 * n41 * n52
                + n15 * n24 * n31 * n43 * n52 - n14 * n25 * n31 * n43 * n52 - n15 * n21 * n34 * n43 * n52
                + n11 * n25 * n34 * n43 * n52 + n14 * n21 * n35 * n43 * n52 - n11 * n24 * n35 * n43 * n52
                - n15 * n23 * n31 * n44 * n52 + n13 * n25 * n31 * n44 * n52 + n15 * n21 * n33 * n44 * n52
                - n11 * n25 * n33 * n44 * n52 - n13 * n21 * n35 * n44 * n52 + n11 * n23 * n35 * n44 * n52
                + n14 * n23 * n31 * n45 * n52 - n13 * n24 * n31 * n45 * n52 - n14 * n21 * n33 * n45 * n52
                + n11 * n24 * n33 * n45 * n52 + n13 * n21 * n34 * n45 * n52 - n11 * n23 * n34 * n45 * n52
                + n15 * n24 * n32 * n41 * n53 - n14 * n25 * n32 * n41 * n53 - n15 * n22 * n34 * n41 * n53
                + n12 * n25 * n34 * n41 * n53 + n14 * n22 * n35 * n41 * n53 - n12 * n24 * n35 * n41 * n53
                - n15 * n24 * n31 * n42 * n53 + n14 * n25 * n31 * n42 * n53 + n15 * n21 * n34 * n42 * n53
                - n11 * n25 * n34 * n42 * n53 - n14 * n21 * n35 * n42 * n53 + n11 * n24 * n35 * n42 * n53
                + n15 * n22 * n31 * n44 * n53 - n12 * n25 * n31 * n44 * n53 - n15 * n21 * n32 * n44 * n53
                + n11 * n25 * n32 * n44 * n53 + n12 * n21 * n35 * n44 * n53 - n11 * n22 * n35 * n44 * n53
                - n14 * n22 * n31 * n45 * n53 + n12 * n24 * n31 * n45 * n53 + n14 * n21 * n32 * n45 * n53
                - n11 * n24 * n32 * n45 * n53 - n12 * n21 * n34 * n45 * n53 + n11 * n22 * n34 * n45 * n53
                - n15 * n23 * n32 * n41 * n54 + n13 * n25 * n32 * n41 * n54 + n15 * n22 * n33 * n41 * n54
                - n12 * n25 * n33 * n41 * n54 - n13 * n22 * n35 * n41 * n54 + n12 * n23 * n35 * n41 * n54
                + n15 * n23 * n31 * n42 * n54 - n13 * n25 * n31 * n42 * n54 - n15 * n21 * n33 * n42 * n54
                + n11 * n25 * n33 * n42 * n54 + n13 * n21 * n35 * n42 * n54 - n11 * n23 * n35 * n42 * n54
                - n15 * n22 * n31 * n43 * n54 + n12 * n25 * n31 * n43 * n54 + n15 * n21 * n32 * n43 * n54
                - n11 * n25 * n32 * n43 * n54 - n12 * n21 * n35 * n43 * n54 + n11 * n22 * n35 * n43 * n54
                + n13 * n22 * n31 * n45 * n54 - n12 * n23 * n31 * n45 * n54 - n13 * n21 * n32 * n45 * n54
                + n11 * n23 * n32 * n45 * n54 + n12 * n21 * n33 * n45 * n54 - n11 * n22 * n33 * n45 * n54
                + n14 * n23 * n32 * n41 * n55 - n13 * n24 * n32 * n41 * n55 - n14 * n22 * n33 * n41 * n55
                + n12 * n24 * n33 * n41 * n55 + n13 * n22 * n34 * n41 * n55 - n12 * n23 * n34 * n41 * n55
                - n14 * n23 * n31 * n42 * n55 + n13 * n24 * n31 * n42 * n55 + n14 * n21 * n33 * n42 * n55
                - n11 * n24 * n33 * n42 * n55 - n13 * n21 * n34 * n42 * n55 + n11 * n23 * n34 * n42 * n55
                + n14 * n22 * n31 * n43 * n55 - n12 * n24 * n31 * n43 * n55 - n14 * n21 * n32 * n43 * n55
                + n11 * n24 * n32 * n43 * n55 + n12 * n21 * n34 * n43 * n55 - n11 * n22 * n34 * n43 * n55
                - n13 * n22 * n31 * n44 * n55 + n12 * n23 * n31 * n44 * n55 + n13 * n21 * n32 * n44 * n55
                - n11 * n23 * n32 * n44 * n55 - n12 * n21 * n33 * n44 * n55 + n11 * n22 * n33 * n44 * n55);

    // ****************************************************
    // **** Loop over six components (beta1-5, fprod)
    // ****************************************************

    fptype temp = smallTerm * (1.0 - Spr0) / (rMassSq - Spr0);

    // Running total
    fpcomplex F0_all(0.0, 0.0);

    fpcomplex _fr11prod(fprod11r, fprod11i);
    fpcomplex _fr12prod(fprod12r, fprod12i);
    fpcomplex _fr13prod(fprod13r, fprod13i);
    fpcomplex _fr14prod(fprod14r, fprod14i);
    fpcomplex _fr15prod(fprod15r, fprod15i);

    fpcomplex _beta_amp;

    pc.incrementIndex(1, 37, 3, 0, 1);

    for(int term = 0; term < 6; term++) { // 0: fprod, 1-5: beta1-5

        // exit loop if this component is not included
        fpcomplex F0(0.0, 0.0); // running sum for each loop

        if(term >= 1) { // beta factors

            if(term == 1)
                _beta_amp = fpcomplex(beta1_ampr, beta1_ampi);
            else if(term == 2)
                _beta_amp = fpcomplex(beta2_ampr, beta2_ampi);
            else if(term == 3)
                _beta_amp = fpcomplex(beta3_ampr, beta3_ampi);
            else if(term == 4)
                _beta_amp = fpcomplex(beta4_ampr, beta4_ampi);
            else if(term == 5)
                _beta_amp = fpcomplex(beta5_ampr, beta5_ampi);

            for(int j = 0; j < 5; j++) {
                fptype tmp = 1. / (poleMassesSq[term - 1] - rMassSq);
                if(fabs(poleMassesSq[term - 1] - rMassSq) < prec) {
                    if(j == 0)
                        F0 += inv11 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 1)
                        F0 += inv12 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 2)
                        F0 += inv13 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 3)
                        F0 += inv14 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 4)
                        F0 += inv15 * g0_matrix[5 * j + (term - 1)];
                } else {
                    tmp *= smallTerm;
                    if(j == 0)
                        F0 += tmp * inv11 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 1)
                        F0 += tmp * inv12 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 2)
                        F0 += tmp * inv13 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 3)
                        F0 += tmp * inv14 * g0_matrix[5 * j + (term - 1)];
                    else if(j == 4)
                        F0 += tmp * inv15 * g0_matrix[5 * j + (term - 1)];
                }
            }
            F0 *= _beta_amp; // don't forget to multiply by amplitude!
        }

        else { // Fprod: changed definition since f11 no longer multiplies all components.
            F0 += inv11 * _fr11prod;
            F0 += inv12 * _fr12prod;
            F0 += inv13 * _fr13prod;
            F0 += inv14 * _fr14prod;
            F0 += inv15 * _fr15prod;
            F0 *= temp;
        }

        F0_all += F0;
    }

    return F0_all;
}

// kMatrixFunction

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
