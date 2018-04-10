/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

This file includes some lineshapes and spinfactors.
Also right now it is the home to some helper functions needed and an implementation of a simple 4-vec class that works
on the GPU
*/

#include <goofit/PDFs/physics/LineshapesPdf.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Version.h>

#if GOOFIT_KMATRIX && THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <goofit/detail/compute_inverse5.h>
#endif

#include <utility>

#include <Eigen/Core>
#include <Eigen/LU>

#include <goofit/detail/Macros.h>

#define NPOLES 5
#define NCHANNELS 5

#define mPiPlus 0.139570
#define mKPlus 0.493677
#define mEta 0.547862
#define mEtap 0.96778

namespace GooFit {

// Lineshape base

// Form factors as in pdg http://pdg.lbl.gov/2012/reviews/rpp2012-rev-dalitz-analysis-formalism.pdf
__device__ fptype BL_PRIME(fptype z2, fptype z02, int L) {
    if(0 == L)
        return 1.0;
    else if(1 == L)
        return (1 + z02) / (1 + z2);
    else if(2 == L)
        return (z02 * z02 + 3 * z02 + 9) / (z2 * z2 + 3 * z2 + 9);
    else {
        printf("ERROR! Oribtal > 2 not supported!\n");
        return 0;
    }

    // Spin 3 and up not accounted for.
}

__device__ fptype BL(fptype z2, int L) {
    if(0 == L)
        return 1.0;
    else if(1 == L)
        return 2 * z2 / (1 + z2);
    else if(2 == L)
        return (13 * z2 * z2) / (z2 * z2 + 3 * z2 + 9);
    else {
        printf("ERROR! Oribtal > 2 not supported!\n");
        return 0;
    }

    // Spin 3 and up not accounted for.
}

__device__ fptype BL2(fptype z2, int L) {
    if(0 == L)
        return 1.0;
    else if(1 == L)
        return 1.0 / (1 + z2);
    else if(2 == L)
        return 1.0 / (z2 * z2 + 3 * z2 + 9);
    else {
        printf("ERROR! Oribtal > 2 not supported!\n");
        return 0;
    }

    // Spin 3 and up not accounted for.
}

__device__ fpcomplex LS_ONE(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) { return {1., 0.}; }

// This function is modeled after BW_BW::getVal() in BW_BW.cpp from the MINT package written by Jonas Rademacker.
__device__ fpcomplex BW(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype meson_radius  = functorConstants[indices[7]];
    fptype resmass       = cudaArray[indices[2]];
    fptype reswidth      = cudaArray[indices[3]];
    unsigned int orbital = indices[4];
    unsigned int FF      = indices[6];

    const unsigned int to2Lplus1 = 2 * orbital + 1;

    fptype mass          = resmass;
    fptype width         = reswidth;
    fptype mumsRecoMass2 = Mpair * Mpair;

    fptype mpsq        = (m1 + m2) * (m1 + m2);
    fptype mmsq        = (m1 - m2) * (m1 - m2);
    fptype num         = (mumsRecoMass2 - mpsq) * (mumsRecoMass2 - mmsq);
    fptype num2        = (mass * mass - mpsq) * (mass * mass - mmsq);
    fptype pABSq       = num / (4 * mumsRecoMass2);
    fptype prSqForGofM = num2 / (4 * mass * mass);
    fptype prSq2       = prSqForGofM < 0 ? 0 : prSqForGofM;
    prSqForGofM        = fabs(prSqForGofM);

    fptype pratio = sqrt(pABSq / prSqForGofM);

    fptype pratio_to_2Jplus1 = 1;

    for(int i = 0; i < to2Lplus1; i++) {
        pratio_to_2Jplus1 *= pratio;
    }

    fptype mratio   = mass / Mpair;
    fptype r2       = meson_radius * meson_radius;
    fptype thisFR   = BL_PRIME(pABSq * r2, prSqForGofM * r2, orbital);
    fptype frFactor = 1;

    if(0 != orbital && 0 != FF) {
        frFactor = (FF == 1 ? BL(pABSq * r2, orbital) : BL_PRIME(pABSq * r2, prSq2 * r2, orbital));
        frFactor = (FF == 3 ? BL2(pABSq * r2, orbital) : frFactor);
    }

    fptype GofM = width * pratio_to_2Jplus1 * mratio * thisFR;

    fptype gamma = mass * sqrt((mass * mass + width * width));
    fptype k     = (2.0 * sqrt(2.0) / M_PI) * mass * width * gamma
               / sqrt(mass * mass + gamma); // Note added additional factor of 2*sqrt(2)/PI here so results are
                                            // comparable to MINT3. MINT2 doesn't have include this.

    fpcomplex BW(mass * mass - mumsRecoMass2, mass * GofM);
    fptype den = (mass * mass - mumsRecoMass2) * (mass * mass - mumsRecoMass2) + mass * GofM * mass * GofM;

    fpcomplex ret = (sqrt(k * frFactor)) / den * BW;
    // printf("m1, m2, Mpair, to2Lplus1, GofM, thisFR, pratio, mratio, pABSq , prSqForGofM, FF, ret.real, ret.imag\n");
    // printf("BW %.7g, %.7g, %.7g, %i, %i, %i, %i\n",meson_radius, resmass, reswidth, orbital, FF, indices[2],
    // indices[3]);
    // printf("BW %.7g, %.7g, %.7g, %i, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g\n", m1, m2, Mpair,
    // to2Lplus1, GofM, thisFR, pratio, mratio, pABSq, prSqForGofM, frFactor, ret.real, ret.imag );
    return ret;
}

// This function is modeled after SBW from the MINT package written by Jonas Rademacker.
__device__ fpcomplex SBW(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype resmass       = GOOFIT_GET_PARAM(2);
    fptype reswidth      = GOOFIT_GET_PARAM(3);
    unsigned int orbital = GOOFIT_GET_INT(4);
    // GOOFIT_GET_INT(5, Mpair, "Mpair");
    unsigned int FF     = GOOFIT_GET_INT(6);
    fptype meson_radius = GOOFIT_GET_CONST(7);

    // fptype meson_radius  = functorConstants[indices[7]];
    // fptype resmass       = cudaArray[indices[2]];
    // fptype reswidth      = cudaArray[indices[3]];
    // unsigned int orbital = indices[4];
    // unsigned int FF      = indices[6];

    fptype mass          = resmass;
    fptype width         = reswidth;
    fptype mumsRecoMass2 = Mpair * Mpair;

    fptype mpsq  = (m1 + m2) * (m1 + m2);
    fptype mmsq  = (m1 - m2) * (m1 - m2);
    fptype num   = (mumsRecoMass2 - mpsq) * (mumsRecoMass2 - mmsq);
    fptype num2  = (mass * mass - mpsq) * (mass * mass - mmsq);
    fptype pABSq = num / (4 * mumsRecoMass2);
    fptype prSq  = num2 / (4 * mass * mass);
    fptype prSq2 = prSq < 0 ? 0 : prSq;
    prSq         = fabs(prSq);

    fptype r2       = meson_radius * meson_radius;
    fptype frFactor = 1;

    if(0 != orbital && 0 != FF) {
        frFactor = (FF == 1 ? BL(pABSq * r2, orbital) : BL_PRIME(pABSq * r2, prSq2 * r2, orbital));
        frFactor = (FF == 3 ? BL2(pABSq * r2, orbital) : frFactor);
    }

    fptype GofM = width;

    fptype gamma = sqrt(mass * mass * (mass * mass + width * width));
    fptype k     = mass * width * gamma / sqrt(mass * mass + gamma);

    fpcomplex BW(mass * mass - mumsRecoMass2, mass * GofM);
    fptype den = (mass * mass - mumsRecoMass2) * (mass * mass - mumsRecoMass2) + mass * GofM * mass * GofM;

    fpcomplex ret = (sqrt(k * frFactor)) / den * BW;

    // printf("m1, m2, Mpair, GofM, pABSq , prSq, FF, ret.real, ret.imag\n");
    // printf("SBW %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g\n", m1, m2, Mpair, GofM, pABSq, prSq, frFactor,
    // ret.real, ret.imag );
    return ret;
}

__device__ fpcomplex bugg_rho2(const fptype &s, const fptype m) {
    fptype rho_squared  = 1. - 4. * m * m / s;
    fpcomplex returnVal = (rho_squared >= 0) ? fpcomplex(1, 0) : fpcomplex(0, 1);
    rho_squared         = (rho_squared >= 0) ? sqrt(rho_squared) : sqrt(-rho_squared);
    return rho_squared * returnVal;
}

__device__ fptype bugg_j1(const fptype &s, const fptype m) {
    fptype rho_pipi  = bugg_rho2(s, m).real();
    fptype returnVal = 2.;
    returnVal += (rho_pipi > 0.) ? rho_pipi * log((1. - rho_pipi) / (1. + rho_pipi)) : 0;
    return returnVal / M_PI;
}

__device__ fptype bugg_Gamma_4pi(const fptype &s,
                                 const fptype mpi,
                                 const fptype &g_4pi,
                                 const fptype &M,
                                 const fptype &lambda_4pi,
                                 const fptype &s0_4pi) {
    fptype returnVal = (s < (16. * mpi * mpi)) ? 0
                                               : g_4pi * (1. / (1 + exp(lambda_4pi * (s0_4pi - s))))
                                                     / (1. / (1 + exp(lambda_4pi * (s0_4pi - M * M))));
    return returnVal;
}

// This function is an adaptation from the bugg lineshape implemented in the MINT package written by Jonas Rademacker.
// this lineshape is not tested yet!
__device__ fpcomplex bugg_MINT(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype s = Mpair * Mpair;

    fptype M          = 0.953;
    fptype b1         = 1.302;
    fptype b2         = 0.340;
    fptype A          = 2.426;
    fptype g_4pi      = 0.011;
    fptype g_2K       = 0.6;
    fptype g_2eta     = 0.2;
    fptype alpha      = 1.3;
    fptype sA         = 0.41;
    fptype s0_4pi     = 7.082 / 2.845;
    fptype lambda_4pi = 2.845;

    fptype g1sq = (b1 + b2 * s) * exp(-(s - M * M) / A);
    fptype z    = bugg_j1(s, mPiPlus) - bugg_j1(M * M, mPiPlus);

    fpcomplex gamma_2pi = fpcomplex(
        g1sq * (s - sA * mPiPlus * mPiPlus) / (M * M - sA * mPiPlus * mPiPlus) * bugg_rho2(s, mPiPlus).real(), 0);
    fpcomplex gamma_2K = g_2K * g1sq * s / (M * M)
                         * exp((-1) * alpha * sqrt((s - 4. * mKPlus * mKPlus) * (s - 4. * mKPlus * mKPlus)))
                         * bugg_rho2(s, mKPlus);
    fpcomplex gamma_2eta = g_2eta * g1sq * s / (M * M)
                           * exp((-1) * alpha * sqrt((s - 4. * mEta * mEta) * (s - 4. * mEta * mEta)))
                           * bugg_rho2(s, mEta);
    fpcomplex gamma_4pi = fpcomplex(bugg_Gamma_4pi(s, mPiPlus, g_4pi, M, lambda_4pi, s0_4pi), 0);

    fpcomplex Gamma_tot = gamma_2pi + gamma_2K + gamma_2eta + gamma_4pi;

    // fpcomplex num = M * gamma_2pi; //only for elastic scattering, not production
    fpcomplex den
        = fpcomplex(M * M - s - M * g1sq * (s - sA * mPiPlus * mPiPlus) / (M * M - sA * mPiPlus * mPiPlus) * z, 0)
          - fpcomplex(0, 1) * M * Gamma_tot;
    fpcomplex returnVal = 1.0 / den;
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",gamma_2pi.real, gamma_2pi.imag, gamma_2K.real,
    // gamma_2K.imag, gamma_2eta.real, gamma_2eta.imag, gamma_4pi.real, gamma_4pi.imag);
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Mpair, Gamma_tot.real, Gamma_tot.imag, g1sq, z,
    // den.real, den.imag, returnVal.real, returnVal.imag);

    // the factor sqrt(1000) gives the correct result in comparison with mint2, I think its because BW/SBW
    // have a factor of sqrt(k) which these lineshapes dont have. For now this stays because it works. further
    // investigation needed.
    return returnVal * sqrt(1000.0);
}

__device__ fpcomplex bugg_MINT3(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype s          = Mpair * Mpair;
    fptype M          = 0.953;
    fptype b1         = 1.302;
    fptype b2         = 0.340;
    fptype A          = 2.426;
    fptype g_4pi      = 0.011;
    fptype g_2K       = 0.6;
    fptype g_2eta     = 0.2;
    fptype alpha      = 1.3;
    fptype s0_4pi     = 7.082 / 2.845;
    fptype lambda_4pi = 2.845;
    fptype sA         = 0.41 * mPiPlus * mPiPlus;

    fptype g1sq      = M * (b1 + b2 * s) * exp(-(s - M * M) / A);
    fptype z         = bugg_j1(s, mPiPlus) - bugg_j1(M * M, mPiPlus);
    fptype adlerZero = (s - sA) / (M * M - sA);

    fptype mk4  = 4. * mKPlus * mKPlus;
    fptype me4  = 4. * mEta * mEta;
    fptype tmp1 = s > mk4 ? s - mk4 : mk4 - s;
    fptype tmp2 = s > me4 ? s - me4 : me4 - s;

    fpcomplex gamma_2pi  = fpcomplex(g1sq * adlerZero * bugg_rho2(s, mPiPlus).real(), 0);
    fpcomplex gamma_2K   = g_2K * g1sq * s / (M * M) * exp((-1) * alpha * tmp1) * bugg_rho2(s, mKPlus);
    fpcomplex gamma_2eta = g_2eta * g1sq * s / (M * M) * exp((-1) * alpha * tmp2) * bugg_rho2(s, mEta);
    fpcomplex gamma_4pi  = fpcomplex(bugg_Gamma_4pi(s, mPiPlus, g_4pi, M, lambda_4pi, s0_4pi), 0);

    fpcomplex Gamma_tot = gamma_2pi + gamma_2K + gamma_2eta + gamma_4pi;

    // fpcomplex num = M * gamma_2pi; //only for elastic scattering, not production
    fpcomplex den       = fpcomplex(M * M - s - adlerZero * g1sq * z, 0) - fpcomplex(0, 1) * Gamma_tot;
    fpcomplex returnVal = 1.0 / den;
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",gamma_2pi.real, gamma_2pi.imag, gamma_2K.real,
    // gamma_2K.imag, gamma_2eta.real, gamma_2eta.imag, gamma_4pi.real, gamma_4pi.imag);
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Mpair, Gamma_tot.real, Gamma_tot.imag, g1sq, z,
    // den.real, den.imag, returnVal.real, returnVal.imag);

    return returnVal;
}

__device__ fpcomplex lass_MINT(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype resmass  = cudaArray[indices[2]];
    fptype reswidth = cudaArray[indices[3]];
    fptype rMass2   = Mpair * Mpair;

    fptype a = 2.07;
    fptype r = 3.32;

    fptype mpsq  = (m1 + m2) * (m1 + m2);
    fptype mmsq  = (m1 - m2) * (m1 - m2);
    fptype num   = (rMass2 - mpsq) * (rMass2 - mmsq);
    fptype num2  = (resmass * resmass - mpsq) * (resmass * resmass - mmsq);
    fptype pABSq = num / (4 * rMass2);
    fptype prSq  = fabs(num2 / (4 * resmass * resmass));

    fptype y          = 2.0 * a * sqrt(pABSq);
    fptype x          = 2.0 + a * r * pABSq;
    fptype cotDeltaBg = x / y;
    fpcomplex phaseshift((cotDeltaBg * cotDeltaBg - 1) / (1 + cotDeltaBg * cotDeltaBg),
                         2 * cotDeltaBg / (1 + cotDeltaBg * cotDeltaBg));
    // (cotDeltaBg*cotDeltaBg-1)/(1+cotDeltaBg*cotDeltaBg) = cos(2*delta)     2*cotDeltaBg / ( 1 +
    // cotDeltaBg*cotDeltaBg) = sin(2*delta)
    fpcomplex den(sqrt(pABSq) * cotDeltaBg, (-1.) * sqrt(pABSq));
    fptype SF           = Mpair * sqrt(prSq) / (resmass * resmass * reswidth);
    fpcomplex BG        = SF / den;
    fpcomplex returnVal = BG + phaseshift * BW(Mpair, m1, m2, indices);
    // printf("Lass: %.5g %.5g %.5g %.5g %.5g %.5g\n",BG.real, BG.imag, phaseshift.real, phaseshift.imag,
    // returnVal.real, returnVal.imag);

    return returnVal;
}

// generalized lass lineshape as implemented in MINT3 by Tim Evans. if F=R=1 and phiF=phiR=0 this is equal to normal
// lass as implemented in Mint3.
// The difference between this and lass mint is not quite clear to me. need to get back to this later.
__device__ fpcomplex glass_MINT3(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype meson_radius  = functorConstants[indices[7]];
    fptype resmass       = cudaArray[indices[2]];
    fptype reswidth      = cudaArray[indices[3]];
    unsigned int orbital = indices[4];
    fptype rMass2        = Mpair * Mpair;

    // fptype a = 2.07;
    // fptype r = 3.32;
    // fptype phiF = 0.0;
    // fptype phiR = 0.0;
    // fptype F = 1.0;
    fptype a    = cudaArray[indices[8]];
    fptype r    = cudaArray[indices[9]];
    fptype phiF = cudaArray[indices[10]];
    fptype phiR = cudaArray[indices[11]];
    fptype F    = cudaArray[indices[12]];

    fptype R = 1.0;
    // printf("GLass: %.5g %.5g %.5g %.5g %.5g %.5g\n",a, r, phiF, phiR, F, R);
    // printf("GLass2: %.5g %.5g %.5g %u \n",meson_radius, resmass, reswidth, orbital);

    fptype mpsq  = (m1 + m2) * (m1 + m2);
    fptype mmsq  = (m1 - m2) * (m1 - m2);
    fptype num   = (rMass2 - mpsq) * (rMass2 - mmsq);
    fptype num2  = (resmass * resmass - mpsq) * (resmass * resmass - mmsq);
    fptype pABSq = num / (4 * rMass2);
    fptype prSq  = fabs(num2 / (4 * resmass * resmass));

    fptype pratio = sqrt(pABSq / prSq);

    fptype pratio_to_2Jplus1 = 1;

    for(int i = 0; i < 2 * orbital + 1; i++) {
        pratio_to_2Jplus1 *= pratio;
    }

    fptype mratio = resmass / Mpair;
    fptype r2     = meson_radius * meson_radius;
    fptype thisFR = BL_PRIME(pABSq * r2, prSq * r2, orbital);
    fptype GofM   = reswidth * pratio_to_2Jplus1 * mratio * thisFR;

    fptype y          = 2.0 * a * sqrt(pABSq);
    fptype x          = 2.0 + a * r * pABSq;
    fptype scattphase = phiF + atan(y / x);
    fptype resphase   = phiR + atan(resmass * GofM / (resmass * resmass - rMass2));
    fptype rho        = 1.0 / sqrt(pABSq / rMass2);
    fpcomplex returnVal
        = (F * sin(scattphase) * fpcomplex(cos(scattphase), sin(scattphase))
           + R * sin(resphase) * fpcomplex(cos(resphase + 2 * scattphase), sin(resphase + 2 * scattphase)))
          * rho;
    // printf("GLass3: %.5g %.5g %.5g %.5g %.5g %.5g\n",rMass2, pABSq, rho, GofM, scattphase, resphase);

    // printf("GLass4: %.5g %.5g\n",returnVal.real, returnVal.imag);
    return returnVal;
}

__device__ fpcomplex aSqrtTerm(const fptype &m0, const fptype &m) {
    fptype a2           = 1 - (2 * m0 / m) * (2 * m0 / m);
    fpcomplex returnVal = a2 > 0 ? fpcomplex(sqrt(a2), 0) : fpcomplex(0, sqrt(-a2));
    return returnVal;
}

__device__ fpcomplex Flatte_MINT(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype meson_radius  = functorConstants[indices[7]];
    fptype resmass       = cudaArray[indices[2]];
    unsigned int orbital = indices[4];
    fptype frFactor      = 1;
    fptype rMass2        = Mpair * Mpair;

    // As far as I understand, this is only valid for the f980
    fptype gPi       = .165;
    fptype gK_by_gPi = 4.21;
    fptype gK        = gPi * gK_by_gPi;
    fptype mPi0      = .1349766;
    fptype mK0       = .497648;

    fptype mpsq = (m1 + m2) * (m1 + m2);
    fptype mmsq = (m1 - m2) * (m1 - m2);
    fptype num  = (rMass2 - mpsq) * (rMass2 - mmsq);
    // fptype num2  = (resmass*resmass - mpsq)*(resmass*resmass - mmsq);
    fptype pABSq = num / (4 * rMass2);
    // fptype prSq = fabs(num2/(4*resmass*resmass));

    fpcomplex Gpipi       = (1. / 3.) * aSqrtTerm(mPi0, Mpair) + (2. / 3.) * aSqrtTerm(mPiPlus, Mpair);
    fpcomplex GKK         = (1. / 2.) * aSqrtTerm(mK0, Mpair) + (1. / 2.) * aSqrtTerm(mKPlus, Mpair);
    fpcomplex FlatteWidth = gPi * Gpipi + gK * GKK;
    // printf("%.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Gpipi.real, Gpipi.imag, GKK.real, GKK.imag, FlatteWidth.real,
    // FlatteWidth.imag, Mpair, pABSq);

    frFactor     = BL2(pABSq * meson_radius * meson_radius, orbital);
    fpcomplex BW = sqrt(frFactor) / fpcomplex(resmass * resmass - rMass2, 0) - fpcomplex(0, 1) * resmass * FlatteWidth;
    return BW;
}

__device__ __thrust_forceinline__ fptype Q2(fptype Msq, fptype M1sq, fptype M2sq) {
    return (Msq / 4. - (M1sq + M2sq) / 2. + (M1sq - M2sq) * (M1sq - M2sq) / (4 * Msq));
}

__device__ __thrust_forceinline__ fptype BlattWeisskopf_Norm(const fptype z2, const fptype z02, unsigned int L) {
    if(L == 0)
        return 1;
    else if(L == 1)
        return (1 + z02) / (1 + z2);
    else if(L == 2)
        return (z02 * z02 + 3 * z02 + 9) / (z2 * z2 + 3 * z2 + 9);
    else {
        abort(__FILE__, __LINE__, "Wrong value of L");
        return 0; // Can't reach
    }
}

__device__ fptype getSpline(fptype x, bool continued, unsigned int *indices) {
    const fptype s_min       = GOOFIT_GET_CONST(8);
    const fptype s_max       = GOOFIT_GET_CONST(9);
    const unsigned int nBins = GOOFIT_GET_INT(10);

    // 11 is the first spine knot, 11+nBins is the first curvature

    if(x <= s_min)
        return continued ? GOOFIT_GET_PARAM(11 + 0) : 0;
    if(x >= s_max)
        return continued ? GOOFIT_GET_PARAM(11 + nBins - 1) : 0;

    fptype spacing = (s_max - s_min) / (nBins - 1.);
    fptype dx      = fmod((x - s_min), spacing);

    auto bin = static_cast<unsigned int>((x - s_min) / spacing);

    fptype m_x_0  = GOOFIT_GET_PARAM(11 + bin);
    fptype m_x_1  = GOOFIT_GET_PARAM(11 + bin + 1);
    fptype m_xf_0 = GOOFIT_GET_CONST(11 + bin + nBins);
    fptype m_xf_1 = GOOFIT_GET_CONST(11 + bin + nBins + 1);

    return m_x_0 + dx * ((m_x_1 - m_x_0) / spacing - (m_xf_1 + 2 * m_xf_0) * spacing / 6) + dx * dx * m_xf_0
           + dx * dx * dx * (m_xf_1 - m_xf_0) / (6 * spacing);
}

__device__ fptype kFactor(fptype mass, fptype width) {
    fptype gamma = mass * sqrt(POW2(mass) + POW2(width));
    fptype k     = 2 * sqrt(2.) * mass * width * gamma / (M_PI * sqrt(POW2(mass) + gamma));
    return sqrt(k);
}

__device__ fpcomplex Spline_TDP(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    const fptype mass  = GOOFIT_GET_PARAM(2);
    const fptype width = GOOFIT_GET_PARAM(3);
    // const unsigned int L = GOOFIT_GET_INT(4);
    const fptype radius = GOOFIT_GET_CONST(7);

    fptype s  = POW2(Mpair);
    fptype s1 = POW2(m1);
    fptype s2 = POW2(m2);

    // This is GSpline.EFF in AmpGen

    fptype q2 = fabs(Q2(s, s1, s2));

    // Non-EFF
    // fptype BF             = sqrt( BlattWeisskopf_Norm(q2 * POW2(radius), 0, L));
    fptype BF = exp(-q2 * POW2(radius) / 2);

    fptype width_shape = width * getSpline(s, true, indices);
    fptype width_norm  = width * getSpline(POW2(mass), false, indices);

    fptype norm          = kFactor(mass, width) * BF;
    fptype running_width = width * width_shape / width_norm;
    fpcomplex iBW        = fpcomplex(POW2(mass) - s, -mass * running_width);
    return norm / iBW;
}

__device__ fpcomplex nonres_DP(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype meson_radius  = functorConstants[indices[7]];
    unsigned int orbital = indices[4];

    fptype mumsRecoMass2 = Mpair * Mpair;

    fptype mpsq       = (m1 + m2) * (m1 + m2);
    fptype mmsq       = (m1 - m2) * (m1 - m2);
    fptype num        = (mumsRecoMass2 - mpsq) * (mumsRecoMass2 - mmsq);
    fptype pABSq      = num / (4 * mumsRecoMass2);
    fptype formfactor = sqrt(BL2(pABSq * meson_radius * meson_radius, orbital));
    // printf("NonRes q2:%.7g FF:%.7g, s %.7g m1 %.7g m2 %.7g r %.7g L %u \n",pABSq, formfactor, mumsRecoMass2,
    // m1,m2,meson_radius, orbital );
    return fpcomplex(1., 0.) * formfactor;
}

__device__ fptype phsp_twoBody(fptype s, fptype m0, fptype m1) { return sqrt(1. - POW2(m0 + m1) / s); }

__device__ fptype phsp_fourPi(fptype s) {
    if(s > 1)
        return phsp_twoBody(s, 2 * mPiPlus, 2 * mPiPlus);
    else
        return 0.00051 + -0.01933 * s + 0.13851 * s * s + -0.20840 * s * s * s + -0.29744 * s * s * s * s
               + 0.13655 * s * s * s * s * s + 1.07885 * s * s * s * s * s * s;
}

#if GOOFIT_KMATRIX

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

__device__ fpcomplex kMatrixFunction(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    // const fptype mass  = GOOFIT_GET_PARAM(2);
    // const fptype width = GOOFIT_GET_PARAM(3);
    // const unsigned int L = GOOFIT_GET_INT(4);
    // const fptype radius = GOOFIT_GET_CONST(7);

    // const fptype pTerm = GOOFIT_GET_INT();

    unsigned int pterm = GOOFIT_GET_INT(8);
    bool is_pole       = GOOFIT_GET_INT(9) == 1;

    fptype sA0      = GOOFIT_GET_PARAM(10);
    fptype sA       = GOOFIT_GET_PARAM(11);
    fptype s0_prod  = GOOFIT_GET_PARAM(12);
    fptype s0_scatt = GOOFIT_GET_PARAM(13);

    Eigen::Array<fptype, NCHANNELS, 1> fscat;
    Eigen::Array<fptype, NPOLES, 1> pmasses;
    Eigen::Array<fptype, NPOLES, NPOLES> couplings;

    for(int i = 0; i < NCHANNELS; i++) {
        fscat(i) = GOOFIT_GET_PARAM(14 + i);
    }

    for(int i = 0; i < NPOLES; i++) {
        for(int j = 0; j < NPOLES; j++)
            couplings(i, j) = GOOFIT_GET_PARAM(14 + NCHANNELS + i * (NPOLES + 1) + j);
        pmasses(i) = GOOFIT_GET_PARAM(14 + NCHANNELS + i * (NPOLES + 1) + NPOLES);
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
#endif

__device__ fptype phsp_FOCUS(fptype s, fptype m0, fptype m1) {
    fptype mp = (m0 + m1);
    fptype mm = (m0 - m1);
    fptype a2 = (1.0 - mp * mp / s) * (1.0 - mm * mm / s);
    return sqrt(a2);
}

__device__ fpcomplex FOCUSFunction(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fptype s         = POW2(Mpair);
    unsigned int mod = GOOFIT_GET_INT(8);

    // mKPlus, mPiPlus, mEtap
    constexpr fptype sNorm = mKPlus * mKPlus + mPiPlus * mPiPlus;
    constexpr fptype s12   = 0.23;
    constexpr fptype s32   = .27;

    fptype I12_adler = (s - s12) / sNorm;
    fptype I32_adler = (s - s32) / sNorm;

    fptype rho1 = phsp_FOCUS(s, mKPlus, mPiPlus);
    fptype rho2 = phsp_FOCUS(s, mKPlus, mEtap);

    fptype pmass = 1.7919;
    Eigen::Array<fptype, 2, 1> coupling;
    coupling << 0.31072, -0.02323;
    fptype X = s / sNorm - 1;

    // constructKMatrix
    Eigen::Array<fptype, 2, 2> kMatrix;

    kMatrix(0, 0) = coupling(0) * coupling(0) / (pmass - s) + 0.79299 - 0.15099 * X + 0.00811 * POW2(X);
    kMatrix(1, 1) = coupling(1) * coupling(1) / (pmass - s) + 0.17054 - 0.0219 * X + 0.00085655 * POW2(X);
    kMatrix(1, 0) = coupling(1) * coupling(0) / (pmass - s) + 0.15040 - 0.038266 * X + 0.0022596 * POW2(X);
    kMatrix(0, 1) = coupling(0) * coupling(1) / (pmass - s) + 0.15040 - 0.038266 * X + 0.0022596 * POW2(X);

    fptype K11 = I12_adler * kMatrix(0, 0);
    fptype K12 = I12_adler * kMatrix(0, 1);
    fptype K22 = I12_adler * kMatrix(1, 1);

    fptype K32 = I32_adler * (-0.22147 + 0.026637 * X - 0.00092057 * POW2(X));

    fptype detK = K11 * K22 - K12 * K12;
    fpcomplex del{1 - rho1 * rho2 * detK, -(rho1 * K11 + rho2 * K22)};

    fpcomplex T11{1., -rho2 * K22};
    fpcomplex T22{1., -rho1 * K11};
    fpcomplex T12{0., rho2 * K12};

    fpcomplex T32 = 1. / fpcomplex(1, -K32 * rho1);

    if(mod == static_cast<unsigned int>(Lineshapes::FOCUS::Mod::Kpi))
        return fpcomplex(K11, -rho2 * detK) / del;
    else if(mod == static_cast<unsigned int>(Lineshapes::FOCUS::Mod::KEta))
        return K12 / del;
    else /*if(mod==Lineshapes::FOCUS::Mod::I32)*/
        return T32;

    return {0., 0.};
}

__device__ resonance_function_ptr ptr_to_LS_ONE     = LS_ONE;
__device__ resonance_function_ptr ptr_to_BW_DP4     = BW;
__device__ resonance_function_ptr ptr_to_lass       = lass_MINT;
__device__ resonance_function_ptr ptr_to_glass3     = glass_MINT3;
__device__ resonance_function_ptr ptr_to_bugg_MINT  = bugg_MINT;
__device__ resonance_function_ptr ptr_to_bugg_MINT3 = bugg_MINT3;
__device__ resonance_function_ptr ptr_to_SBW        = SBW;
__device__ resonance_function_ptr ptr_to_NONRES_DP  = nonres_DP;
__device__ resonance_function_ptr ptr_to_Flatte     = Flatte_MINT;
__device__ resonance_function_ptr ptr_to_Spline     = Spline_TDP;
__device__ resonance_function_ptr ptr_to_FOCUS      = FOCUSFunction;
#if GOOFIT_KMATRIX
__device__ resonance_function_ptr ptr_to_kMatrix = kMatrixFunction;
#endif

// This constructor is protected
Lineshape::Lineshape(std::string name)
    : GooPdf(name) {
    // Making room for index of decay-related constants. Assumption:
    // These are mother mass and three daughter masses in that order.
    // They will be registered by the object that uses this resonance,
    // which will tell this object where to find them by calling setConstantIndex.
}

std::vector<fptype> make_spline_curvatures(std::vector<Variable> vars, Lineshapes::spline_t SplineInfo) {
    size_t size = std::get<2>(SplineInfo) - 2;
    Eigen::Matrix<fptype, Eigen::Dynamic, Eigen::Dynamic> m(size, size);
    for(size_t i = 0; i < size; i++) {
        m(i, i) = 4;
        if(i != size - 1) {
            m(i, i + 1) = 1;
            m(i + 1, 1) = 1;
        }
    }
    m = m.inverse();

    Eigen::Matrix<fptype, Eigen::Dynamic, 1> L(size);
    for(unsigned int i = 0; i < size; ++i)
        L[i] = vars[i + 2].getValue() - 2 * vars[i + 1].getValue() + vars[i].getValue();

    auto mtv       = m * L;
    fptype spacing = (std::get<0>(SplineInfo) - std::get<1>(SplineInfo)) / std::get<2>(SplineInfo);

    std::vector<fptype> ret(vars.size(), 0);
    for(unsigned int i = 0; i < size; ++i) {
        ret.at(i + 1) = 6 * mtv(i) / POW2(spacing);
    }
    return ret;
}

Lineshapes::GSpline::GSpline(std::string name,
                             Variable mass,
                             Variable width,
                             unsigned int L,
                             unsigned int Mpair,
                             FF FormFac,
                             fptype radius,
                             std::vector<Variable> AdditionalVars,
                             spline_t SplineInfo)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    if(std::get<2>(SplineInfo) != AdditionalVars.size())
        throw GeneralError("bins {} != vars {}", std::get<2>(SplineInfo), AdditionalVars.size());
    GOOFIT_ADD_CONST(8, std::get<0>(SplineInfo), "MinSpline");
    GOOFIT_ADD_CONST(9, std::get<1>(SplineInfo), "MaxSpline");
    GOOFIT_ADD_INT(10, std::get<2>(SplineInfo), "NSpline");

    int i = 11;
    for(auto &par : AdditionalVars) {
        GOOFIT_ADD_PARAM(i++, par, "Knot");
    }

    std::vector<fptype> SplineCTerms = make_spline_curvatures(AdditionalVars, SplineInfo);

    for(auto &par : SplineCTerms) {
        // Calc curve
        GOOFIT_ADD_CONST(i++, par, "CKnot");
    }

    GET_FUNCTION_ADDR(ptr_to_Spline);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::GLASS::GLASS(std::string name,
                         Variable mass,
                         Variable width,
                         unsigned int L,
                         unsigned int Mpair,
                         FF FormFac,
                         fptype radius,
                         std::vector<Variable> AdditionalVars)
    : Lineshape(name) {
    if(5 != AdditionalVars.size()) {
        throw GeneralError("It seems you forgot to provide the vector with the five necessary variables for GLASS, a, "
                           "r, phiF, phiR and F (in that order)");
    }

    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    for(int i = 0; i < 5; i++) {
        GOOFIT_ADD_PARAM(8 + i, AdditionalVars[i], "LassVars");
    }

    GET_FUNCTION_ADDR(ptr_to_glass3);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::One::One(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_LS_ONE);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::RBW::RBW(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_BW_DP4);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::LASS::LASS(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_lass);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::NonRes::NonRes(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_NONRES_DP);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::Bugg::Bugg(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_bugg_MINT);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::Bugg3::Bugg3(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_bugg_MINT3);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::Flatte::Flatte(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_Flatte);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::SBW::SBW(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GET_FUNCTION_ADDR(ptr_to_SBW);

    GOOFIT_FINALIZE_PDF;
}

Lineshapes::FOCUS::FOCUS(std::string name,
                         Mod mod,
                         Variable mass,
                         Variable width,
                         unsigned int L,
                         unsigned int Mpair,
                         FF FormFac,
                         fptype radius)
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GOOFIT_ADD_INT(8, static_cast<unsigned int>(mod), "Lineshape modifier");

    GET_FUNCTION_ADDR(ptr_to_FOCUS);

    GOOFIT_FINALIZE_PDF;
}

#if GOOFIT_KMATRIX
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
    : Lineshape(name) {
    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    GOOFIT_ADD_INT(8, pterm, "pTerm");
    GOOFIT_ADD_INT(9, is_pole ? 1 : 0, "Channel");

    GOOFIT_ADD_PARAM(10, sA0, "sA0");
    GOOFIT_ADD_PARAM(11, sA, "sA");
    GOOFIT_ADD_PARAM(12, s0_prod, "s0_prod");
    GOOFIT_ADD_PARAM(13, s0_scatt, "s0_scatt");

    for(int i = 0; i < NCHANNELS; i++) {
        GOOFIT_ADD_PARAM(14 + i, fscat.at(i), "fscat");
    }

    for(int i = 0; i < NPOLES * (NPOLES + 1); i++) {
        GOOFIT_ADD_PARAM(14 + NCHANNELS + i, poles.at(i), "Poles");
    }

    GET_FUNCTION_ADDR(ptr_to_kMatrix);

    GOOFIT_FINALIZE_PDF;
}
#endif

Amplitude::Amplitude(std::string uniqueDecayStr,
                     Variable ar,
                     Variable ai,
                     std::vector<Lineshape *> LS,
                     std::vector<SpinFactor *> SF,
                     unsigned int nPerm)
    : _uniqueDecayStr(std::move(uniqueDecayStr))
    , _ar(ar)
    , _ai(ai)
    , _SF(std::move(SF))
    , _LS(std::move(LS))
    , _nPerm(nPerm) {}

bool Amplitude::operator==(const Amplitude &A) const {
    return _uniqueDecayStr == A._uniqueDecayStr && _ar == A._ar && _ai == A._ai && _LS == A._LS && _SF == A._SF
           and _nPerm == A._nPerm;
}

} // namespace GooFit
