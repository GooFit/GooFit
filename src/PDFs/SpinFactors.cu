/*
04/05/2016 Christoph Hasse
The Spinfactors are an adaptation from the MINT implementation, by Jonas Rademacker.

DISCLAIMER:
This code is not sufficently tested yet and still under heavy development!
*/
#include "goofit/PDFs/physics/SpinFactors.h"
#include "goofit/PDFs/physics/SpinHelper.h"

namespace GooFit {

#define ZEMACH 1

__device__ fptype DtoPP1_PtoSP2_StoP3P4(fptype *Vecs, unsigned int *indices) { return 1.0; }

__device__ fptype ONE(fptype *Vecs, unsigned int *indices) { return 1.0; }

__device__ fptype FF_12_34_L1(fptype *Vecs, unsigned int *indices) {
    fptype mother_radius = functorConstants[indices[1]];
    unsigned int p1      = indices[2];
    unsigned int p2      = indices[3];
    unsigned int p3      = indices[4];
    unsigned int p4      = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    fptype m1 = (P1 + P2).Mag2();
    fptype m2 = (P3 + P4).Mag2();
    fptype s  = (P1 + P2 + P3 + P4).Mag2();
    fptype q2 = s / 4. - (m1 + m2) / 2. + (m1 - m2) * (m1 - m2) / (4 * s);
    fptype z2 = q2 * mother_radius * mother_radius;
    fptype ff = 1.0 / (1 + z2);
    // printf("%.5g, %.5g, %.5g, %.5g\n",s,m1,m2,sqrt(ff) );
    return sqrt(ff);
}

__device__ fptype FF_12_34_L2(fptype *Vecs, unsigned int *indices) {
    fptype mother_radius = functorConstants[indices[1]];
    unsigned int p1      = indices[2];
    unsigned int p2      = indices[3];
    unsigned int p3      = indices[4];
    unsigned int p4      = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    fptype m1 = (P1 + P2).Mag2();
    fptype m2 = (P3 + P4).Mag2();
    fptype s  = (P1 + P2 + P3 + P4).Mag2();
    fptype q2 = s / 4. - (m1 + m2) / 2. + (m1 - m2) * (m1 - m2) / (4 * s);
    fptype z2 = q2 * mother_radius * mother_radius;
    fptype ff = 1.0 / (z2 * z2 + 3 * z2 + 9);

    return sqrt(ff);
}

__device__ fptype FF_123_4_L1(fptype *Vecs, unsigned int *indices) {
    fptype mother_radius = functorConstants[indices[1]];
    unsigned int p1      = indices[2];
    unsigned int p2      = indices[3];
    unsigned int p3      = indices[4];
    unsigned int p4      = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    fptype m1 = (P1 + P2 + P3).Mag2();
    fptype m2 = P4.Mag2();
    fptype s  = (P1 + P2 + P3 + P4).Mag2();
    fptype q2 = s / 4. - (m1 + m2) / 2. + (m1 - m2) * (m1 - m2) / (4 * s);
    fptype z2 = q2 * mother_radius * mother_radius;
    fptype ff = 1.0 / (1 + z2);

    return sqrt(ff);
}

__device__ fptype FF_123_4_L2(fptype *Vecs, unsigned int *indices) {
    fptype mother_radius = functorConstants[indices[1]];
    unsigned int p1      = indices[2];
    unsigned int p2      = indices[3];
    unsigned int p3      = indices[4];
    unsigned int p4      = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    fptype m1 = (P1 + P2 + P3).Mag2();
    fptype m2 = P4.Mag2();
    fptype s  = (P1 + P2 + P3 + P4).Mag2();
    fptype q2 = s / 4. - (m1 + m2) / 2. + (m1 - m2) * (m1 - m2) / (4 * s);
    fptype z2 = q2 * mother_radius * mother_radius;
    fptype ff = 1.0 / (z2 * z2 + 3 * z2 + 9);

    return sqrt(ff);
}

__device__ fptype DtoPP1_PtoVP2_VtoP3P4(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pV = P3 + P4;
    gpuLVec qV = P3 - P4;

#ifdef ZEMACH
    gpuLVec pP = pV + P2;
    gpuLVec qP = pV - P2;

    ZTspin1 LP(qP, pP, pP.M());
    ZTspin1 LV(qV, pV, pV.M());

    return LP.Dot(LV);
#else

    return (P2.Dot(qV) - P2.Dot(pV) * pV.Dot(qV)) / (pV.M() * pV.M());
#endif
}

__device__ fptype DtoV1V2_V1toP1P2_V2toP3P4_S(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",0, P1.getX(), P1.getY(), P1.getZ(),P1.getE());
    // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",1, P2.getX(), P2.getY(), P2.getZ(),P2.getE());
    // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",2, P3.getX(), P3.getY(), P3.getZ(),P3.getE());
    // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",3, P4.getX(), P4.getY(), P4.getZ(),P4.getE());

    gpuLVec pV1 = P1 + P2;
    gpuLVec qV1 = P1 - P2;
    gpuLVec pV2 = P3 + P4;
    gpuLVec qV2 = P3 - P4;

    fptype MV1 = sqrt(pV1.Dot(pV1));
    fptype MV2 = sqrt(pV2.Dot(pV2));

    fptype returnVal
        = (qV1.Dot(qV2) - qV1.Dot(pV1) * pV1.Dot(qV2) / (MV1 * MV1) - qV1.Dot(pV2) * pV2.Dot(qV2) / (MV2 * MV2)
           + qV1.Dot(pV1) * pV1.Dot(pV2) * pV2.Dot(qV2) / (MV1 * MV1 * MV2 * MV2));
    // printf("s1 %.5g; %i,%i,%i,%i\n",returnVal, indices[2], indices[3], indices[4], indices[5]);
    return returnVal;
}

__device__ fptype DtoV1V2_V1toP1P2_V2toP3P4_P(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pV1 = P1 + P2;
    gpuLVec qV1 = P2 - P1;
    gpuLVec pV2 = P3 + P4;
    gpuLVec qV2 = P4 - P3;

    gpuLVec pD = pV1 + pV2;
    gpuLVec qD = pV1 - pV2;

#ifdef ZEMACH
    fptype MV1 = sqrt(pV1.Dot(pV1));
    fptype MV2 = sqrt(pV2.Dot(pV2));

    ZTspin1 LD(qD, pD, pD.M());
    ZTspin1 LV1(qV1, pV1, MV1);
    ZTspin1 LV2(qV2, pV2, MV2);

    return -LeviCivita(pD, LD, LV1).Dot(LV2); // minus gives the same result as MINT3

#else
    return LeviCivita(pD, qD, qV1, qV2);
#endif
}

__device__ fptype DtoV1V2_V1toP1P2_V2toP3P4_D(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pV1 = P1 + P2;
    gpuLVec qV1 = P1 - P2;
    gpuLVec pV2 = P3 + P4;
    gpuLVec qV2 = P3 - P4;

#ifdef ZEMACH
    gpuLVec pD = pV1 + pV2;
    gpuLVec qD = pV1 - pV2;
    double mD  = pD.M();

    ZTspin1 tV1(qV1, pV1, pV1.M());
    ZTspin1 tV2(qV2, pV2, pV2.M());
    ZTspin2 tD(qD, pD, mD);
    double returnVal = tV1.Contract(tD.Contract(tV2));
    // printf("%f, %f, %f, %f,",P1.GetX(), P1.GetY(), P1.GetZ(), P1.GetE() );
    // printf("%f, %f, %f, %f,",P2.GetX(), P2.GetY(), P2.GetZ(), P2.GetE() );
    // printf("%f, %f, %f, %f,",P3.GetX(), P3.GetY(), P3.GetZ(), P3.GetE() );
    // printf("%f, %f, %f, %f,",P4.GetX(), P4.GetY(), P4.GetZ(), P4.GetE() );
    // printf("%f\n",returnVal );
    return returnVal;
#else

    fptype MV1       = sqrt(pV1.Dot(pV1));
    fptype MV2       = sqrt(pV2.Dot(pV2));
    fptype returnVal = (qV1.Dot(pV2) - qV1.Dot(pV1) * pV1.Dot(pV2) / (MV1 * MV1))
                       * (qV2.Dot(pV1) - qV2.Dot(pV2) * pV2.Dot(pV1) / (MV2 * MV2));

    return returnVal;
#endif
}

__device__ fptype DtoV1P1_V1toV2P2_V2toP3P4(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pV1 = P2 + P3 + P4;
    gpuLVec pV2 = P3 + P4;
    gpuLVec qV1 = (P3 + P4) - P2;
    gpuLVec qV2 = P3 - P4;

#ifdef ZEMACH
    double MV1 = pV1.M();
    double MV2 = pV2.M();

    gpuLVec pD = pV1 + P1;
    gpuLVec qD = pV1 - P1;

    ZTspin1 LD(qD, pD, pD.M());
    ZTspin1 LV1(qV1, pV1, MV1);
    ZTspin1 LV2(qV2, pV2, MV2);
    SpinSumV PV1(pV1, MV1);

    gpuLVec tmp = PV1.Dot(LD);

    return LeviCivita(pV1, LV1, LV2).Dot(tmp);
#else

    // printf("%f, %f, %f, %f\n",P1.getX(), P1.getY(), P1.getZ(), P1.getE() );
    // printf("%f, %f, %f, %f\n",P2.getX(), P2.getY(), P2.getZ(), P2.getE() );
    // printf("%f, %f, %f, %f\n",P3.getX(), P3.getY(), P3.getZ(), P3.getE() );
    // printf("%f, %f, %f, %f\n",P4.getX(), P4.getY(), P4.getZ(), P4.getE() );

    fptype returnVal = LeviCivita(pV1, qV1, P1, qV2);
    return returnVal;
#endif
}

__device__ fptype DtoVS_VtoP1P2_StoP3P4(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pS = P3 + P4;
    gpuLVec pV = P1 + P2;
    gpuLVec qV = P1 - P2;

#ifdef ZEMACH
    gpuLVec pD = pV + pS;
    gpuLVec qD = pV - pS;
    ZTspin1 LD(qD, pD, pD.M());
    ZTspin1 LV(qV, pV, pV.M());
    return (LD.Dot(LV));
#else

    // printf("%f, %f, %f, %f\n",P1.getX(), P1.getY(), P1.getZ(), P1.getE() );
    // printf("%f, %f, %f, %f\n",P2.getX(), P2.getY(), P2.getZ(), P2.getE() );
    // printf("%f, %f, %f, %f\n",P3.getX(), P3.getY(), P3.getZ(), P3.getE() );
    // printf("%f, %f, %f, %f\n",P4.getX(), P4.getY(), P4.getZ(), P4.getE() );

    fptype MV = sqrt(pV.Dot(pV));

    fptype returnVal = (pS.Dot(qV) - pS.Dot(pV) * pV.Dot(qV) / (MV * MV));
    return returnVal;
#endif
}

__device__ fptype DtoAP1_AtoSP2_StoP3P4(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pS = P3 + P4;
    gpuLVec pA = P2 + pS;
    gpuLVec qA = pS - P2;
    gpuLVec pD = pA + P1;
    gpuLVec qD = pA - P1;

// printf("%f, %f, %f, %f\n",P1.GetX(), P1.GetY(), P1.GetZ(), P1.GetE() );
// printf("%f, %f, %f, %f\n",P2.GetX(), P2.GetY(), P2.GetZ(), P2.GetE() );
// printf("%f, %f, %f, %f\n",P3.GetX(), P3.GetY(), P3.GetZ(), P3.GetE() );
// printf("%f, %f, %f, %f\n",P4.GetX(), P4.GetY(), P4.GetZ(), P4.GetE() );

#ifdef ZEMACH
    ZTspin1 LD(qD, pD, pD.M());
    ZTspin1 LA(qA, pA, pA.M());
    return (LD.Dot(LA));
#else

    fptype MA        = sqrt(pA.Dot(pA));
    fptype returnVal = (P1.Dot(qA) - P1.Dot(pA) * pA.Dot(qA) / (MA * MA));
    return returnVal;
#endif
}

__device__ fptype DtoAP1_AtoVP2_VtoP3P4(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pV = P3 + P4;
    gpuLVec qV = P3 - P4;
    gpuLVec pA = P2 + pV;
    gpuLVec pD = P1 + pA;
    gpuLVec qD = pA - P1;

#ifdef ZEMACH
    ZTspin1 LB(qD, pD, pD.M());
    ZTspin1 LV(qV, pV, pV.M());
    SpinSumV PA(pA, pA.M());

    gpuLVec tmp = PA.Dot(LV);
    return -(LB.Dot(tmp)); // minus to be equal to MINT3
#else
    gpuLVec p0 = P1;

    fptype MA        = sqrt(pA.Dot(pA));
    fptype MV        = sqrt(pV.Dot(pV));
    fptype returnVal = P1.Dot(qV) - p0.Dot(pA) * pA.Dot(qV) / (MA * MA) - p0.Dot(pV) * pV.Dot(qV) / (MV * MV)
                       + p0.Dot(pA) * pA.Dot(pV) * pV.Dot(qV) / (MA * MA * MV * MV);
    // printf("spin %.7g\n",returnVal );
    return returnVal;
#endif
}

__device__ fptype DtoAP1_AtoVP2Dwave_VtoP3P4(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pV = P3 + P4;
    gpuLVec qV = P3 - P4;
    gpuLVec pA = P2 + pV;
    gpuLVec qA = P2 - pV;
    gpuLVec pD = P1 + pA;
    gpuLVec qD = pA - P1;

#ifdef ZEMACH
    ZTspin1 LD(qD, pD, pD.M());
    ZTspin1 LV(qV, pV, pV.M());
    ZTspin2 LA(qA, pA, pA.M());
    gpuLVec tmp = LA.Contract(LV);

    return (LD.Dot(tmp));
#else

    fptype MA = sqrt(pA.Dot(pA));
    fptype MV = sqrt(pV.Dot(pV));

    fptype returnVal
        = P1.Dot((qA - qA.Dot(pA) * pA * (1. / (MA * MA)))) * (qV - qV.Dot(pV) * pV * (1. / (MV * MV))).Dot(pA);

    // printf("spin %.7g\n",returnVal );
    return returnVal;
#endif
}

__device__ fptype DtoTP1_TtoVP2_VtoP3P4(fptype *Vecs, unsigned int *indices) {
    unsigned int p1 = indices[2];
    unsigned int p2 = indices[3];
    unsigned int p3 = indices[4];
    unsigned int p4 = indices[5];
    gpuLVec P1(Vecs[0 + 4 * p1], Vecs[1 + 4 * p1], Vecs[2 + 4 * p1], Vecs[3 + 4 * p1]);
    gpuLVec P2(Vecs[0 + 4 * p2], Vecs[1 + 4 * p2], Vecs[2 + 4 * p2], Vecs[3 + 4 * p2]);
    gpuLVec P3(Vecs[0 + 4 * p3], Vecs[1 + 4 * p3], Vecs[2 + 4 * p3], Vecs[3 + 4 * p3]);
    gpuLVec P4(Vecs[0 + 4 * p4], Vecs[1 + 4 * p4], Vecs[2 + 4 * p4], Vecs[3 + 4 * p4]);

    gpuLVec pV = P3 + P4;
    gpuLVec qV = P3 - P4;
    gpuLVec pT = P2 + pV;
    gpuLVec qT = pV - P2;
    gpuLVec pD = P1 + pT;
    gpuLVec qD = pT - P1;

    ZTspin2 t2T(qT, pT, pT.M());
    ZTspin1 tV(qV, pV, pV.M());

// gpuLVec DT(t2T.Contract(qD));

#ifdef ZEMACH
    ZTspin1 tD(qD, pD, pD.M());
    ZTspin1 t1T(qT, pT, pT.M());

    SpinSumV P1T(pT, pT.M());
    SpinSumT P2T(pT, pT.M());

    gpuLVec tmp = LeviCivita(t1T, pT, P1T.Dot(tV));
    return P2T.Sandwich(tD, tD, tmp, t1T) + 1. / 3. * tD.Dot(tD) / (pT.Dot(pT)) * P2T.Sandwich(pD, pD, tmp, t1T);
#else

    fptype returnVal = LeviCivita(pD, qD, DT, tV);

    // printf("spin %.7g\n",returnVal );
    return returnVal;
#endif
}

__device__ spin_function_ptr ptr_to_DtoPP1_PtoSP2_StoP3P4       = DtoPP1_PtoSP2_StoP3P4;
__device__ spin_function_ptr ptr_to_DtoPP1_PtoVP2_VtoP3P4       = DtoPP1_PtoVP2_VtoP3P4;
__device__ spin_function_ptr ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_S = DtoV1V2_V1toP1P2_V2toP3P4_S;
__device__ spin_function_ptr ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_P = DtoV1V2_V1toP1P2_V2toP3P4_P;
__device__ spin_function_ptr ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_D = DtoV1V2_V1toP1P2_V2toP3P4_D;
__device__ spin_function_ptr ptr_to_DtoVS_VtoP1P2_StoP3P4       = DtoVS_VtoP1P2_StoP3P4;
__device__ spin_function_ptr ptr_to_DtoV1P1_V1toV2P2_V2toP3P4   = DtoV1P1_V1toV2P2_V2toP3P4;
__device__ spin_function_ptr ptr_to_DtoAP1_AtoSP2_StoP3P4       = DtoAP1_AtoSP2_StoP3P4;
__device__ spin_function_ptr ptr_to_DtoAP1_AtoVP2_VtoP3P4       = DtoAP1_AtoVP2_VtoP3P4;
__device__ spin_function_ptr ptr_to_DtoAP1_AtoVP2Dwave_VtoP3P4  = DtoAP1_AtoVP2Dwave_VtoP3P4;
__device__ spin_function_ptr ptr_to_DtoTP1_TtoVP2_VtoP3P4       = DtoTP1_TtoVP2_VtoP3P4;

__device__ spin_function_ptr ptr_to_FF_12_34_L1 = FF_12_34_L1;
__device__ spin_function_ptr ptr_to_FF_12_34_L2 = FF_12_34_L2;
__device__ spin_function_ptr ptr_to_FF_123_4_L1 = FF_123_4_L1;
__device__ spin_function_ptr ptr_to_FF_123_4_L2 = FF_123_4_L2;
__device__ spin_function_ptr ptr_to_ONE         = ONE;

SpinFactor::SpinFactor(
    std::string name, SF_4Body SF, unsigned int P0, unsigned int P1, unsigned int P2, unsigned int P3)
    : GooPdf(nullptr, name)
    , _SF(SF)
    , _P0(P0)
    , _P1(P1)
    , _P2(P2)
    , _P3(P3) {
    std::vector<unsigned int> pindices;
    pindices.push_back(0); // dummy for index to constants.
    pindices.push_back(P0);
    pindices.push_back(P1);
    pindices.push_back(P2);
    pindices.push_back(P3);

    switch(SF) {
    case SF_4Body::DtoPP1_PtoSP2_StoP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoPP1_PtoSP2_StoP3P4);
        break;

    case SF_4Body::DtoPP1_PtoVP2_VtoP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoPP1_PtoVP2_VtoP3P4);
        break;

    case SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S:
        GET_FUNCTION_ADDR(ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_S);
        break;

    case SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P:
        GET_FUNCTION_ADDR(ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_P);
        break;

    case SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D:
        GET_FUNCTION_ADDR(ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_D);
        break;

    case SF_4Body::DtoAP1_AtoVP2_VtoP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoAP1_AtoVP2_VtoP3P4);
        break;

    case SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoAP1_AtoVP2Dwave_VtoP3P4);
        break;

    case SF_4Body::DtoVS_VtoP1P2_StoP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoVS_VtoP1P2_StoP3P4);
        break;

    case SF_4Body::DtoV1P1_V1toV2P2_V2toP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoV1P1_V1toV2P2_V2toP3P4);
        break;

    case SF_4Body::DtoAP1_AtoSP2_StoP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoAP1_AtoSP2_StoP3P4);
        break;

    case SF_4Body::DtoTP1_TtoVP2_VtoP3P4:
        GET_FUNCTION_ADDR(ptr_to_DtoTP1_TtoVP2_VtoP3P4);
        break;

    case SF_4Body::FF_12_34_L1:
        GET_FUNCTION_ADDR(ptr_to_FF_12_34_L1);
        break;

    case SF_4Body::FF_12_34_L2:
        GET_FUNCTION_ADDR(ptr_to_FF_12_34_L2);
        break;

    case SF_4Body::FF_123_4_L1:
        GET_FUNCTION_ADDR(ptr_to_FF_123_4_L1);
        break;

    case SF_4Body::FF_123_4_L2:
        GET_FUNCTION_ADDR(ptr_to_FF_123_4_L2);
        break;

    case SF_4Body::ONE:
        GET_FUNCTION_ADDR(ptr_to_ONE);
        break;

    default:
        std::cout << "No Spinfunction implemented for that kind." << std::endl;
        exit(0);
        break;
    }

    initialize(pindices);
}

} // namespace GooFit
