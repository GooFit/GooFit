/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

This file includes some lineshapes and spinfactors.
Also right now it is the home to some helper functions needed and an implementation of a simple 4-vec class that works
on the GPU
*/

#include "goofit/PDFs/physics/SpinFactors.h"
#include "goofit/PDFs/physics/LineshapesPdf.h"

#include <utility>

#include <goofit/detail/Macros.h>

namespace GooFit {

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

__device__ fpcomplex LS_ONE(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) { return fpcomplex(1, 0); }

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

    if(0 != orbital and 0 != FF) {
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

    if(0 != orbital and 0 != FF) {
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
    fptype mPiPlus    = .13957018;
    fptype mKPlus     = .493677;
    fptype mEta       = .54751;

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
    fptype mPiPlus    = .13957018;
    fptype mKPlus     = .493677;
    fptype mEta       = .547862;
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
    fptype mPiPlus   = .13957018;
    fptype mKPlus    = .493677;
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

__device__ fpcomplex Spline_TDP(fptype Mpair, fptype m1, fptype m2, unsigned int *indices) {
    fpcomplex ret(0, 0);

    const unsigned int nKnobs = indices[4]; // orbital

    unsigned int idx                 = 5; // Next index
    unsigned int i                   = 0;
    const unsigned int pwa_coefs_idx = idx;
    idx += 2 * nKnobs;
    const fptype *mKKlimits = &(functorConstants[indices[idx]]);
    fptype mAB = Mpair, mAC = m1, mBC = m2;

    int khiAB = 0, khiAC = 0;
    fptype dmKK, aa, bb, aa3, bb3;
    unsigned int timestorun = 1; //+doSwap;
    while(khiAB < nKnobs) {
        if(mAB < mKKlimits[khiAB])
            break;
        khiAB++;
    }

    if(khiAB <= 0 || khiAB == nKnobs)
        timestorun = 0;
    while(khiAC < nKnobs) {
        if(mAC < mKKlimits[khiAC])
            break;
        khiAC++;
    }

    if(khiAC <= 0 || khiAC == nKnobs)
        timestorun = 0;

    for(i = 0; i < timestorun; i++) {
        unsigned int kloAB                = khiAB - 1; //, kloAC = khiAC -1;
        unsigned int twokloAB             = kloAB + kloAB;
        unsigned int twokhiAB             = khiAB + khiAB;
        fptype pwa_coefs_real_kloAB       = cudaArray[indices[pwa_coefs_idx + twokloAB]];
        fptype pwa_coefs_real_khiAB       = cudaArray[indices[pwa_coefs_idx + twokhiAB]];
        fptype pwa_coefs_imag_kloAB       = cudaArray[indices[pwa_coefs_idx + twokloAB + 1]];
        fptype pwa_coefs_imag_khiAB       = cudaArray[indices[pwa_coefs_idx + twokhiAB + 1]];
        fptype pwa_coefs_prime_real_kloAB = 0; // cDeriatives[twokloAB];
        fptype pwa_coefs_prime_real_khiAB = 0; // cDeriatives[twokhiAB];
        fptype pwa_coefs_prime_imag_kloAB = 0; // cDeriatives[twokloAB+1];
        fptype pwa_coefs_prime_imag_khiAB = 0; // cDeriatives[twokhiAB+1];
        //  printf("m12: %f: %f %f %f %f %f %f %d %d %d\n", mAB, mKKlimits[0], mKKlimits[nKnobs-1],
        //  pwa_coefs_real_khiAB, pwa_coefs_imag_khiAB, pwa_coefs_prime_real_khiAB, pwa_coefs_prime_imag_khiAB, khiAB,
        //  khiAC, timestorun );

        dmKK = mKKlimits[khiAB] - mKKlimits[kloAB];
        aa   = (mKKlimits[khiAB] - mAB) / dmKK;
        bb   = 1 - aa;
        aa3  = aa * aa * aa;
        bb3  = bb * bb * bb;
        //  ret += aa * pwa_coefs[kloAB] + bb * pwa_coefs[khiAB] + ((aa3 - aa)*pwa_coefs_prime[kloAB] + (bb3 - bb) *
        //  pwa_coefs_prime[khiAB]) * (dmKK*dmKK)/6.0;
        ret.real(ret.real() + aa * pwa_coefs_real_kloAB + bb * pwa_coefs_real_khiAB
                 + ((aa3 - aa) * pwa_coefs_prime_real_kloAB + (bb3 - bb) * pwa_coefs_prime_real_khiAB) * (dmKK * dmKK)
                       / 6.0);
        ret.imag(ret.imag() + aa * pwa_coefs_imag_kloAB + bb * pwa_coefs_imag_khiAB
                 + ((aa3 - aa) * pwa_coefs_prime_imag_kloAB + (bb3 - bb) * pwa_coefs_prime_imag_khiAB) * (dmKK * dmKK)
                       / 6.0);
        khiAB = khiAC;
        mAB   = mAC;
    }
    return ret;
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
    return fpcomplex(1, 0) * formfactor;
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

Lineshape::Lineshape(std::string name,
                     Variable *mass,
                     Variable *width,
                     unsigned int L,
                     unsigned int Mpair,
                     LS kind,
                     FF FormFac,
                     fptype radius,
                     std::vector<Variable *> AdditionalVars,
                     spline_t SplineInfo)
    : GooPdf(nullptr, name)
    , _mass(mass)
    , _width(width)
    , _L(L)
    , _Mpair(Mpair)
    , _kind(kind)
    , _FormFac(FormFac)
    , _SplineInfo(SplineInfo) {
    GOOFIT_START_PDF;

    // Making room for index of decay-related constants. Assumption:
    // These are mother mass and three daughter masses in that order.
    // They will be registered by the object that uses this resonance,
    // which will tell this object where to find them by calling setConstantIndex.

    GOOFIT_ADD_PARAM(2, mass, "mass");
    GOOFIT_ADD_PARAM(3, width, "width");

    GOOFIT_ADD_INT(4, L, "L");
    GOOFIT_ADD_INT(5, Mpair, "Mpair");

    GOOFIT_ADD_INT(6, enum_to_underlying(FormFac), "FormFac");

    GOOFIT_ADD_CONST(7, radius, "radius");

    // pindices.push_back(registerParameter(mass));
    // pindices.push_back(registerParameter(width));
    // pindices.push_back(L);
    // pindices.push_back(Mpair);
    // pindices.push_back(enum_to_underlying(FormFac));

    // pindices.push_back(registerConstants(1));
    // MEMCPY_TO_SYMBOL(functorConstants, &radius, sizeof(fptype), cIndex * sizeof(fptype), cudaMemcpyHostToDevice);

    switch(kind) {
    case LS::ONE:
        GET_FUNCTION_ADDR(ptr_to_LS_ONE);
        break;

    case LS::BW:
        GET_FUNCTION_ADDR(ptr_to_BW_DP4);
        break;

    case LS::Lass:
        GET_FUNCTION_ADDR(ptr_to_lass);
        break;

    case LS::Lass_M3:
        if(5 != AdditionalVars.size()) {
            fprintf(stderr,
                    "It seems you forgot to provide the vector with the five necessary variables for GLASS, a, "
                    "r, phiF, phiR and F (in that order)");
            exit(0);
        }

        for(int i = 0; i < 5; i++) {
            GOOFIT_ADD_PARAM(8 + i, AdditionalVars[i], "LassVars");
            // pindices.push_back(registerParameter(AdditionalVars[i]));
        }

        GET_FUNCTION_ADDR(ptr_to_glass3);
        break;

    case LS::nonRes:
        GET_FUNCTION_ADDR(ptr_to_NONRES_DP);
        break;

    case LS::Bugg:
        GET_FUNCTION_ADDR(ptr_to_bugg_MINT);
        break;

    case LS::Bugg3:
        GET_FUNCTION_ADDR(ptr_to_bugg_MINT3);
        break;

    case LS::SBW:
        GET_FUNCTION_ADDR(ptr_to_SBW);
        break;

    case LS::Flatte:
        GET_FUNCTION_ADDR(ptr_to_Flatte);
        break;

    case LS::Spline:
        if(std::get<2>(_SplineInfo) != AdditionalVars.size())
            throw GeneralError("bins {} != vars {}", std::get<2>(_SplineInfo), AdditionalVars.size());
        GOOFIT_ADD_CONST(8, std::get<0>(_SplineInfo), "MinSpline");
        GOOFIT_ADD_CONST(9, std::get<1>(_SplineInfo), "MaxSpline");
        GOOFIT_ADD_CONST(10, std::get<2>(_SplineInfo), "NSpline");
        {
            int i = 11;
            for(auto &par : AdditionalVars) {
                GOOFIT_ADD_PARAM(i++, par, "Knot");
            }
        }

        GET_FUNCTION_ADDR(ptr_to_Spline);
        break;

    default:
        throw GeneralError("It seems that the requested lineshape is not implemented yet. Check LineshapesPdf.cu");
    }

    GOOFIT_FINALIZE_PDF;
}

Lineshape::Lineshape(std::string name)
    : GooPdf(nullptr, name) {
    GOOFIT_START_PDF;
    GET_FUNCTION_ADDR(ptr_to_NONRES_DP);
    GOOFIT_FINALIZE_PDF;
}

Amplitude::Amplitude(std::string uniqueDecayStr,
                     Variable *ar,
                     Variable *ai,
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
    return _uniqueDecayStr == A._uniqueDecayStr and _ar == A._ar and _ai == A._ai and _LS == A._LS and _SF == A._SF
           and _nPerm == A._nPerm;
}

} // namespace GooFit
