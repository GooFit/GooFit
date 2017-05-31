/*
04/05/2016 Christoph Hasse
DISCLAIMER:
This code is not sufficently tested yet and still under heavy development!

Helper functions
*/

#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

__device__ fptype Mass(const fptype *P0) {
    return sqrt(-P0[0] * P0[0] - P0[1] * P0[1] - P0[2] * P0[2] + P0[3] * P0[3]);
}
__device__ fptype Mass(const fptype *P0, const fptype *P1) {
    return sqrt(-((P0[0] + P1[0]) * (P0[0] + P1[0])) - ((P0[1] + P1[1]) * (P0[1] + P1[1]))
                - ((P0[2] + P1[2]) * (P0[2] + P1[2]))
                + ((P0[3] + P1[3]) * (P0[3] + P1[3])));
}
__device__ fptype Mass(const fptype *P0, const fptype *P1, const fptype *P2) {
    return sqrt(-((P0[0] + P1[0] + P2[0]) * (P0[0] + P1[0] + P2[0]))
                - ((P0[1] + P1[1] + P2[1]) * (P0[1] + P1[1] + P2[1]))
                - ((P0[2] + P1[2] + P2[2]) * (P0[2] + P1[2] + P2[2]))
                + ((P0[3] + P1[3] + P2[3]) * (P0[3] + P1[3] + P2[3])));
}
__device__ fptype VecDot(const fptype *P0, const fptype *P1) {
    return (P0[0] * P1[0] + P0[1] + P1[1] + P0[2] + P1[2] + P0[3] + P1[3]);
}

__device__ void get4Vecs(fptype *Vecs,
                         const unsigned int &constants,
                         const fptype &m12,
                         const fptype &m34,
                         const fptype &cos12,
                         const fptype &cos34,
                         const fptype &phi) {
    fptype M  = functorConstants[constants + 1];
    fptype m1 = functorConstants[constants + 2];
    fptype m2 = functorConstants[constants + 3];
    fptype m3 = functorConstants[constants + 4];
    fptype m4 = functorConstants[constants + 5];
    // printf("g4v %f, %f, %f, %f, %f\n",M, m1, m2, m3, m4 );
    fptype E1     = (m12 * m12 + m1 * m1 - m2 * m2) / (2 * m12);
    fptype E2     = (m12 * m12 - m1 * m1 + m2 * m2) / (2 * m12);
    fptype E3     = (m34 * m34 + m3 * m3 - m4 * m4) / (2 * m34);
    fptype E4     = (m34 * m34 - m3 * m3 + m4 * m4) / (2 * m34);
    fptype p1     = sqrt(E1 * E1 - m1 * m1);
    fptype p3     = sqrt(E3 * E3 - m3 * m3);
    fptype sin12  = sqrt(1 - cos12 * cos12);
    fptype sin34  = sqrt(1 - cos34 * cos34);
    fptype ED1    = (M * M + m12 * m12 - m34 * m34) / (2 * m12);
    fptype PD1    = sqrt(ED1 * ED1 - M * M);
    fptype beta1  = PD1 / ED1;
    fptype gamma1 = 1.0 / sqrt(1 - beta1 * beta1);
    fptype ED2    = (M * M - m12 * m12 + m34 * m34) / (2 * m34);
    fptype PD2    = sqrt(ED2 * ED2 - M * M);
    fptype beta2  = -PD2 / ED2;
    fptype gamma2 = 1.0 / sqrt(1 - beta2 * beta2);
    // printf("g4v %f, %f, %f, %f, %f\n",E1, m1, E2, p1, p3 );

    // set X-component
    Vecs[0]  = cos12 * p1;
    Vecs[4]  = -cos12 * p1;
    Vecs[8]  = -cos34 * p3;
    Vecs[12] = cos34 * p3;

    // set Y-component
    Vecs[1]  = sin12 * p1;
    Vecs[5]  = -sin12 * p1;
    Vecs[9]  = -sin34 * p3;
    Vecs[13] = sin34 * p3;

    // set Z-component
    Vecs[2]  = 0;
    Vecs[6]  = 0;
    Vecs[10] = 0;
    Vecs[14] = 0;

    // set E-component
    Vecs[3]  = E1;
    Vecs[7]  = E2;
    Vecs[11] = E3;
    Vecs[15] = E4;

    fptype tmpE = Vecs[3];
    fptype tmpX = Vecs[0];
    Vecs[3]     = gamma1 * (tmpE + beta1 * tmpX);
    Vecs[0]     = gamma1 * (tmpX + beta1 * tmpE);

    tmpE    = Vecs[7];
    tmpX    = Vecs[4];
    Vecs[7] = gamma1 * (tmpE + beta1 * tmpX);
    Vecs[4] = gamma1 * (tmpX + beta1 * tmpE);

    tmpE     = Vecs[11];
    tmpX     = Vecs[8];
    Vecs[11] = gamma2 * (tmpE + beta2 * tmpX);
    Vecs[8]  = gamma2 * (tmpX + beta2 * tmpE);

    tmpE     = Vecs[15];
    tmpX     = Vecs[12];
    Vecs[15] = gamma2 * (tmpE + beta2 * tmpX);
    Vecs[12] = gamma2 * (tmpX + beta2 * tmpE);

    // rotation around X-axis of the first two vectors.
    fptype cosphi = cos(phi);
    fptype sinphi = sin(phi);

    // note that Z-component is zero thus rotation is as easy as this:
    Vecs[2] = sinphi * Vecs[1];
    Vecs[1] = cosphi * Vecs[1];

    Vecs[6] = sinphi * Vecs[5];
    Vecs[5] = cosphi * Vecs[5];
}

__device__ fptype getmass(const unsigned int &pair,
                          fptype &d1,
                          fptype &d2,
                          const fptype *vecs,
                          const fptype &m1,
                          const fptype &m2,
                          const fptype &m3,
                          const fptype &m4) {
    const fptype *P1 = vecs;
    const fptype *P2 = (vecs + 4);
    const fptype *P3 = (vecs + 8);
    const fptype *P4 = (vecs + 12);
    fptype mpair     = 0;

    switch(pair) {
    case 2:
        d1    = m1;
        d2    = m3;
        mpair = Mass(P1, P3);
        break;

    case 3:
        d1    = m1;
        d2    = m4;
        mpair = Mass(P1, P4);
        break;

    case 4:
        d1    = m2;
        d2    = m3;
        mpair = Mass(P2, P3);
        break;

    case 5:
        d1    = m2;
        d2    = m4;
        mpair = Mass(P2, P4);
        break;

    case 6:
        d1    = Mass(P1, P2);
        d2    = m3;
        mpair = Mass(P1, P2, P3);
        break;

    case 7:
        d1    = Mass(P1, P3);
        d2    = m2;
        mpair = Mass(P1, P2, P3);
        break;

    case 8:
        d1    = Mass(P2, P3);
        d2    = m1;
        mpair = Mass(P1, P2, P3);
        break;

    case 9:
        d1    = Mass(P1, P2);
        d2    = m4;
        mpair = Mass(P1, P2, P4);
        break;

    case 10:
        d1    = Mass(P1, P4);
        d2    = m2;
        mpair = Mass(P1, P2, P4);
        break;

    case 11:
        d1    = Mass(P2, P4);
        d2    = m1;
        mpair = Mass(P1, P2, P4);
        break;

    case 12:
        d1    = Mass(P1, P3);
        d2    = m4;
        mpair = Mass(P1, P3, P4);
        break;

    case 13:
        d1    = Mass(P1, P4);
        d2    = m3;
        mpair = Mass(P1, P3, P4);
        break;

    case 14:
        d1    = Mass(P3, P4);
        d2    = m1;
        mpair = Mass(P1, P3, P4);
        break;

    case 15:
        d1    = Mass(P2, P3);
        d2    = m4;
        mpair = Mass(P2, P3, P4);
        break;

    case 16:
        d1    = Mass(P2, P4);
        d2    = m3;
        mpair = Mass(P2, P3, P4);
        break;

    case 17:
        d1    = Mass(P3, P4);
        d2    = m2;
        mpair = Mass(P2, P3, P4);
        break;
    }

    return mpair;
}
} // namespace GooFit
