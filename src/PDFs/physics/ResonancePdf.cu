#include "goofit/PDFs/physics/ResonancePdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

//__device__ fpcomplex cDerivaties[1000];

__device__ fptype twoBodyCMmom(double rMassSq, fptype d1m, fptype d2m) {
    // For A -> B + C, calculate momentum of B and C in rest frame of A.
    // PDG 38.16.
    fptype kin1 = 1 - POW2(d1m + d2m) / rMassSq;

    if(kin1 >= 0)
        kin1 = sqrt(kin1);
    else
        kin1 = 1;

    fptype kin2 = 1 - POW2(d1m - d2m) / rMassSq;

    if(kin2 >= 0)
        kin2 = sqrt(kin2);
    else
        kin2 = 1;

    return 0.5 * sqrt(rMassSq) * kin1 * kin2;
}

__device__ fptype dampingFactorSquare(const fptype &cmmom, const int &spin, const fptype &mRadius) {
    fptype square = mRadius * mRadius * cmmom * cmmom;
    fptype dfsq   = 1 + square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    fptype dfsqres = dfsq + 8 + 2 * square + square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}

__device__ fptype spinFactor(unsigned int spin,
                             fptype motherMass,
                             fptype daug1Mass,
                             fptype daug2Mass,
                             fptype daug3Mass,
                             fptype m12,
                             fptype m13,
                             fptype m23,
                             unsigned int cyclic_index) {
    if(0 == spin)
        return 1; // Should not cause branching since every thread evaluates the same resonance at the same time.

    /*
    // Copied from BdkDMixDalitzAmp

    fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
    fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug3Mass));
    fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));

    fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m12 : m12));
    fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m23 : m13));
    fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

    // The above, collapsed into single tests where possible.
    fptype _mA = (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass);
    fptype _mB = (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass);
    fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));

    fptype _mAC = (PAIR_23 == cyclic_index ? m13 : m23);
    fptype _mBC = (PAIR_12 == cyclic_index ? m13 : m12);
    fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    */

    // Copied from EvtDalitzReso, with assumption that pairAng convention matches pipipi0 from EvtD0mixDalitz.
    // Again, all threads should get the same branch.
    fptype _mA  = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
    fptype _mB  = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
    fptype _mC  = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
    fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12));
    fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m12 : m13));
    fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

    fptype massFactor = 1.0 / _mAB;
    fptype sFactor    = -1;
    sFactor *= ((_mBC - _mAC) + (massFactor * (motherMass * motherMass - _mC * _mC) * (_mA * _mA - _mB * _mB)));

    if(2 == spin) {
        sFactor *= sFactor;
        fptype extraterm = ((_mAB - (2 * motherMass * motherMass) - (2 * _mC * _mC))
                            + massFactor * POW2(motherMass * motherMass - _mC * _mC));
        extraterm *= ((_mAB - (2 * _mA * _mA) - (2 * _mB * _mB)) + massFactor * POW2(_mA * _mA - _mB * _mB));
        extraterm /= 3;
        sFactor -= extraterm;
    }

    return sFactor;
}

__device__ fpcomplex plainBW(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    // fptype motherMass   = c_motherMass;//RO_CACHE(pc.constants[pc.constantIdx + 1]);
    // fptype daug1Mass    = c_daug1Mass;//RO_CACHE(pc.constants[pc.constantIdx + 2]);
    // fptype daug2Mass    = c_daug2Mass;//RO_CACHE(pc.constants[pc.constantIdx + 3]);
    // fptype daug3Mass    = c_daug3Mass;//RO_CACHE(pc.constants[pc.constantIdx + 4]);
    // fptype meson_radius = c_meson_radius;//RO_CACHE(pc.constants[pc.constantIdx + 5]);

    unsigned int spin         = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    unsigned int cyclic_index = RO_CACHE(pc.constants[pc.constantIdx + 2]);

    fptype resmass  = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype reswidth = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);

    fptype rMassSq  = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype frFactor = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).
    fptype measureDaughterMoms = twoBodyCMmom(rMassSq,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(resmass,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass));

    if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, c_meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
    }

    // RBW evaluation
    fptype A = (resmass - rMassSq);
    fptype B = resmass * reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor
               / sqrt(rMassSq);
    fptype C = 1.0 / (A * A + B * B);
    fpcomplex ret(A * C, B * C); // Dropping F_D=1

    ret *= sqrt(frFactor);
    fptype spinF = spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);
    ret *= spinF;
    pc.incrementIndex(1, 2, 2, 0, 1);
    // printf("%f, %f, %f, %f\n",ret.real, ret.imag, m12, m13);
    return ret;
}

__device__ fpcomplex gaussian(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    // indices[1] is unused constant index, for consistency with other function types.
    fptype resmass  = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype reswidth = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);

    unsigned int cyclic_index = RO_CACHE(pc.constants[pc.constantIdx + 1]);

    // Notice sqrt - this function uses mass, not mass-squared like the other resonance types.
    fptype massToUse = sqrt(PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    massToUse -= resmass;
    massToUse /= reswidth;
    massToUse *= massToUse;
    fptype ret = exp(-0.5 * massToUse);

    pc.incrementIndex(1, 2, 1, 0, 1);
    // Ignore factor 1/sqrt(2pi).
    ret /= reswidth;

    return fpcomplex(ret, 0);
}

__device__ fptype hFun(double s, double daug2Mass, double daug3Mass) {
    // Last helper function
    const fptype _pi = 3.14159265359;
    double sm        = daug2Mass + daug3Mass;
    double SQRTs     = sqrt(s);
    double k_s       = twoBodyCMmom(s, daug2Mass, daug3Mass);

    double val = ((2 / _pi) * (k_s / SQRTs) * log((SQRTs + 2 * k_s) / (sm)));

    return val;
}

__device__ fptype dh_dsFun(double s, double daug2Mass, double daug3Mass) {
    // Yet another helper function
    const fptype _pi = 3.14159265359;
    double k_s       = twoBodyCMmom(s, daug2Mass, daug3Mass);

    double val = hFun(s, daug2Mass, daug3Mass) * (1.0 / (8.0 * POW2(k_s)) - 1.0 / (2.0 * s)) + 1.0 / (2.0 * _pi * s);
    return val;
}

__device__ fptype dFun(double s, double daug2Mass, double daug3Mass) {
    // Helper function used in Gronau-Sakurai
    const fptype _pi = 3.14159265359;
    double sm        = daug2Mass + daug3Mass;
    double sm24      = sm * sm / 4.0;
    double m         = sqrt(s);
    double k_m2      = twoBodyCMmom(s, daug2Mass, daug3Mass);

    double val = 3.0 / _pi * sm24 / POW2(k_m2) * log((m + 2 * k_m2) / sm) + m / (2 * _pi * k_m2)
                 - sm24 * m / (_pi * POW3(k_m2));
    return val;
}

__device__ fptype fsFun(double s, double m2, double gam, double daug2Mass, double daug3Mass) {
    // Another G-S helper function

    double k_s   = twoBodyCMmom(s, daug2Mass, daug3Mass);
    double k_Am2 = twoBodyCMmom(m2, daug2Mass, daug3Mass);

    double f = gam * m2 / POW3(k_Am2);
    f *= (POW2(k_s) * (hFun(s, daug2Mass, daug3Mass) - hFun(m2, daug2Mass, daug3Mass))
          + (m2 - s) * POW2(k_Am2) * dh_dsFun(m2, daug2Mass, daug3Mass));

    return f;
}

__device__ fpcomplex gouSak(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    // fptype motherMass   = c_motherMass;//RO_CACHE(pc.constants[pc.constantIdx + 1]);
    // fptype daug1Mass    = c_daug1Mass;//RO_CACHE(pc.constants[pc.constantIdx + 2]);
    // fptype daug2Mass    = c_daug2Mass;//RO_CACHE(pc.constants[pc.constantIdx + 3]);
    // fptype daug3Mass    = c_daug3Mass;//RO_CACHE(pc.constants[pc.constantIdx + 4]);
    // fptype meson_radius = c_meson_radius;//RO_CACHE(pc.constants[pc.constantIdx + 5]);

    unsigned int spin         = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    unsigned int cyclic_index = RO_CACHE(pc.constants[pc.constantIdx + 2]);

    fptype resmass  = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype reswidth = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);

    fptype rMassSq  = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype frFactor = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).
    fptype measureDaughterMoms = twoBodyCMmom(rMassSq,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(resmass,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass));

    if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, c_meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
    }

    // Implement Gou-Sak:

    fptype D = (1.0 + dFun(resmass, c_daug2Mass, c_daug3Mass) * reswidth / sqrt(resmass));
    fptype E = resmass - rMassSq + fsFun(rMassSq, resmass, reswidth, c_daug2Mass, c_daug3Mass);
    fptype F = sqrt(resmass) * reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor;

    D /= (E * E + F * F);
    fpcomplex retur(D * E, D * F); // Dropping F_D=1
    retur *= sqrt(frFactor);
    retur *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

    pc.incrementIndex(1, 2, 2, 0, 1);

    return retur;
}

__device__ fpcomplex lass(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    // fptype motherMass   = c_motherMass;//RO_CACHE(pc.constants[pc.constantIdx + 1]);
    // fptype daug1Mass    = c_daug1Mass;//RO_CACHE(pc.constants[pc.constantIdx + 2]);
    // fptype daug2Mass    = c_daug2Mass;//RO_CACHE(pc.constants[pc.constantIdx + 3]);
    // fptype daug3Mass    = c_daug3Mass;//RO_CACHE(pc.constants[pc.constantIdx + 4]);
    // fptype meson_radius = c_meson_radius;//RO_CACHE(pc.constants[pc.constantIdx + 5]);

    unsigned int spin         = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    unsigned int cyclic_index = RO_CACHE(pc.constants[pc.constantIdx + 2]);

    fptype resmass  = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype reswidth = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);

    fptype rMassSq  = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype frFactor = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).

    fptype measureDaughterMoms = twoBodyCMmom(rMassSq,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_23 == cyclic_index ? c_daug3Mass : c_daug2Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(resmass,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_23 == cyclic_index ? c_daug3Mass : c_daug2Mass));

    if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, c_meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
    }

    // Implement LASS:
    /*
    fptype s = kinematics(m12, m13, _trackinfo[i]);
    fptype q = twoBodyCMmom(s, _trackinfo[i]);
    fptype m0  = _massRes[i]->getValFast();
    fptype _g0 = _gammaRes[i]->getValFast();
    int spin   = _spinRes[i];
    fptype g = runningWidthFast(s, m0, _g0, spin, _trackinfo[i], FrEval(s, m0, _trackinfo[i], spin));
    */

    fptype q = measureDaughterMoms;
    fptype g = reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor / sqrt(rMassSq);

    fptype _a    = 0.22357;
    fptype _r    = -15.042;
    fptype _R    = 1; // ?
    fptype _phiR = 1.10644;
    fptype _B    = 0.614463;
    fptype _phiB = -0.0981907;

    // background phase motion
    fptype cot_deltaB  = (1.0 / (_a * q)) + 0.5 * _r * q;
    fptype qcot_deltaB = (1.0 / _a) + 0.5 * _r * q * q;

    // calculate resonant part
    fpcomplex expi2deltaB = fpcomplex(qcot_deltaB, q) / fpcomplex(qcot_deltaB, -q);
    fpcomplex resT        = fpcomplex(cos(_phiR + 2 * _phiB), sin(_phiR + 2 * _phiB)) * _R;

    fpcomplex prop = fpcomplex(1, 0) / fpcomplex(resmass - rMassSq, sqrt(resmass) * g);
    // resT *= prop*m0*_g0*m0/twoBodyCMmom(m0*m0, _trackinfo[i])*expi2deltaB;
    resT *= prop * (resmass * reswidth / nominalDaughterMoms) * expi2deltaB;

    // calculate bkg part
    resT += fpcomplex(cos(_phiB), sin(_phiB)) * _B * (cos(_phiB) + cot_deltaB * sin(_phiB)) * sqrt(rMassSq)
            / fpcomplex(qcot_deltaB, -q);

    resT *= sqrt(frFactor);
    resT *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

    pc.incrementIndex(1, 2, 2, 0, 1);

    return resT;
}

__device__ fpcomplex nonres(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    pc.incrementIndex(1, 0, 0, 0, 1);
    return fpcomplex(1, 0);
}

__device__ void
getAmplitudeCoefficients(fpcomplex a1, fpcomplex a2, fptype &a1sq, fptype &a2sq, fptype &a1a2real, fptype &a1a2imag) {
    // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
    a1sq = thrust::norm(a1);
    a2sq = thrust::norm(a2);
    a1 *= conj(a2);
    a1a2real = a1.real();
    a1a2imag = a1.imag();
}

__device__ fpcomplex flatte(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    // indices[1] is unused constant index, for consistency with other function types.
    fptype resmass            = pc.parameters[pc.parameterIdx + 1];
    fptype g1                 = pc.parameters[pc.parameterIdx + 2];
    fptype g2                 = pc.parameters[pc.parameterIdx + 3] * g1;
    unsigned int cyclic_index = pc.constants[pc.constantIdx + 1];
    unsigned int doSwap       = pc.constants[pc.constantIdx + 2];

    fptype pipmass = 0.13957018;
    fptype pi0mass = 0.1349766;
    fptype kpmass  = 0.493677;
    fptype k0mass  = 0.497614;

    fptype twopimasssq  = 4 * pipmass * pipmass;
    fptype twopi0masssq = 4 * pi0mass * pi0mass;
    fptype twokmasssq   = 4 * kpmass * kpmass;
    fptype twok0masssq  = 4 * k0mass * k0mass;

    fpcomplex ret(0., 0.);
    for(int i = 0; i < 1 + doSwap; i++) {
        fptype rhopipi_real = 0, rhopipi_imag = 0;
        fptype rhokk_real = 0, rhokk_imag = 0;

        fptype s = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

        if(s >= twopimasssq)
            rhopipi_real += (2. / 3) * sqrt(1 - twopimasssq / s); // Above pi+pi- threshold
        else
            rhopipi_imag += (2. / 3) * sqrt(-1 + twopimasssq / s);
        if(s >= twopi0masssq)
            rhopipi_real += (1. / 3) * sqrt(1 - twopi0masssq / s); // Above pi0pi0 threshold
        else
            rhopipi_imag += (1. / 3) * sqrt(-1 + twopi0masssq / s);
        if(s >= twokmasssq)
            rhokk_real += 0.5 * sqrt(1 - twokmasssq / s); // Above K+K- threshold
        else
            rhokk_imag += 0.5 * sqrt(-1 + twokmasssq / s);
        if(s >= twok0masssq)
            rhokk_real += 0.5 * sqrt(1 - twok0masssq / s); // Above K0K0 threshold
        else
            rhokk_imag += 0.5 * sqrt(-1 + twok0masssq / s);
        fptype A = (resmass * resmass - s) + resmass * (rhopipi_imag * g1 + rhokk_imag * g2);
        fptype B = resmass * (rhopipi_real * g1 + rhokk_real * g2);
        fptype C = 1.0 / (A * A + B * B);
        fpcomplex retur(A * C, B * C);
        ret += retur;
        if(doSwap) {
            fptype swpmass = m12;
            m12            = m13;
            m13            = swpmass;
        }
    }

    pc.incrementIndex (1, 3, 2, 0, 1);

    return ret;
}

__device__ fpcomplex cubicspline(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    fpcomplex ret(0, 0);
    unsigned int cyclic_index        = pc.constants[pc.constantIdx + 2];
    unsigned int doSwap              = pc.constants[pc.constantIdx + 3];
    const unsigned int nKnobs        = pc.constants[pc.constantIdx + 4];
    unsigned int idx                 = 5; // Next index
    unsigned int i                   = 0;
    const unsigned int pwa_coefs_idx = idx;
    idx += 2 * nKnobs;
    //const fptype *mKKlimits = &(functorConstants[indices[idx]]);
    fptype mAB = m12, mAC = m13, mBC = m23;
    switch(cyclic_index) {
    case PAIR_13:
        mAB = m13;
        mAC = m12;
        break;
    case PAIR_23:
        mAB = m23;
        mAC = m12;
        mBC = m13;
        break;
    }

    int khiAB = 0, khiAC = 0;
    fptype dmKK, aa, bb, aa3, bb3;
    unsigned int timestorun = 1 + doSwap;
    while(khiAB < nKnobs) {
        if(mAB < pc.constants[pc.constantIdx + khiAB])
            break;
        khiAB++;
    }

    if(khiAB <= 0 || khiAB == nKnobs)
        timestorun = 0;
    while(khiAC < nKnobs) {
        if(mAC < pc.constants[pc.constantIdx + khiAC])
            break;
        khiAC++;
    }

    if(khiAC <= 0 || khiAC == nKnobs)
        timestorun = 0;

    for(i = 0; i < timestorun; i++) {
        unsigned int kloAB                = khiAB - 1; //, kloAC = khiAC -1;
        unsigned int twokloAB             = kloAB + kloAB;
        unsigned int twokhiAB             = khiAB + khiAB;
        fptype pwa_coefs_real_kloAB       = pc.parameters[pc.parameterIdx + pwa_coefs_idx + twokloAB];
        fptype pwa_coefs_real_khiAB       = pc.parameters[pc.parameterIdx + pwa_coefs_idx + twokhiAB];
        fptype pwa_coefs_imag_kloAB       = pc.parameters[pc.parameterIdx + pwa_coefs_idx + twokloAB + 1];
        fptype pwa_coefs_imag_khiAB       = pc.parameters[pc.parameterIdx + pwa_coefs_idx + twokhiAB + 1];
        //fptype pwa_coefs_prime_real_kloAB = cDeriatives[twokloAB];
        fptype pwa_coefs_prime_real_kloAB = 1.0;
        fptype pwa_coefs_prime_real_khiAB = 1.0;
        fptype pwa_coefs_prime_imag_kloAB = 1.0;
        fptype pwa_coefs_prime_imag_khiAB = 1.0;
        //fptype pwa_coefs_prime_real_khiAB = cDeriatives[twokhiAB];
        //fptype pwa_coefs_prime_imag_kloAB = cDeriatives[twokloAB + 1];
        //fptype pwa_coefs_prime_imag_khiAB = cDeriatives[twokhiAB + 1];
        //  printf("m12: %f: %f %f %f %f %f %f %d %d %d\n", mAB, mKKlimits[0], mKKlimits[nKnobs-1],
        //  pwa_coefs_real_khiAB, pwa_coefs_imag_khiAB, pwa_coefs_prime_real_khiAB, pwa_coefs_prime_imag_khiAB, khiAB,
        //  khiAC, timestorun );

        dmKK = pc.constants[pc.constantIdx + khiAB] - pc.constants[pc.constantIdx + kloAB];
        aa   = (pc.constants[pc.constantIdx + khiAB] - mAB) / dmKK;
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

__device__ resonance_function_ptr ptr_to_RBW      = plainBW;
__device__ resonance_function_ptr ptr_to_GOUSAK   = gouSak;
__device__ resonance_function_ptr ptr_to_GAUSSIAN = gaussian;
__device__ resonance_function_ptr ptr_to_NONRES   = nonres;
__device__ resonance_function_ptr ptr_to_LASS     = lass;
__device__ resonance_function_ptr ptr_to_FLATTE   = flatte;
__device__ resonance_function_ptr ptr_to_SPLINE   = cubicspline;

namespace Resonances {

RBW::RBW(std::string name,
        Variable ar,
        Variable ai,
        Variable mass,
        Variable width,
        unsigned int sp,
        unsigned int cyc,
        bool symmDP) : ResonancePdf(name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);
    registerConstant(symmDP);

    resonanceType = 0;
}
 

GS::GS(std::string name,
       Variable ar,
       Variable ai,
       Variable mass,
       Variable width,
       unsigned int sp,
       unsigned int cyc,
       bool symmDP)
    : ResonancePdf(name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);
    registerConstant(symmDP);

    resonanceType = 1;
}

LASS::LASS(std::string name,
           Variable ar,
           Variable ai,
           Variable mass,
           Variable width,
           unsigned int sp,
           unsigned int cyc,
           bool symmDP)
    : ResonancePdf(name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);
    registerConstant(symmDP);

    resonanceType = 2;
}

// Constructor for regular BW,Gounaris-Sakurai,LASS
Gauss::Gauss(std::string name, Variable ar, Variable ai, Variable mass, Variable width, unsigned int cyc)
    : ResonancePdf(name, ar, ai) {
    // Making room for index of decay-related constants. Assumption:
    // These are mother mass and three daughter masses in that order.
    // They will be registered by the object that uses this resonance,
    // which will tell this object where to find them by calling setConstantIndex.
    registerParameter(mass);
    registerParameter(width);

    registerConstant(cyc);

    resonanceType = 3;
}

NonRes::NonRes(std::string name, Variable ar, Variable ai)
    : ResonancePdf(name, ar, ai) {
    resonanceType = 4;
}

FLATTE::FLATTE(std::string name,
               Variable ar,
               Variable ai,
               Variable mean,
               Variable g1,
               Variable rg2og1,
               unsigned int cyc,
               bool symmDP)
    : ResonancePdf(name, ar, ai) {
    registerParameter(mean);
    registerParameter(g1);
    registerParameter(rg2og1);

    registerConstant(cyc);
    registerConstant(symmDP);

    resonanceType = 5;
}

Spline::Spline(std::string name,
               Variable ar,
               Variable ai,
               std::vector<fptype> &HH_bin_limits,
               std::vector<Variable> &pwa_coefs_reals,
               std::vector<Variable> &pwa_coefs_imags,
               unsigned int cyc,
               const bool symmDP)
    : ResonancePdf(name, ar, ai) {
    const unsigned int nKnobs = HH_bin_limits.size();

    registerConstant(cyc);
    registerConstant(symmDP);
    registerConstant(nKnobs);

    for(int i = 0; i < pwa_coefs_reals.size(); i++) {
        registerConstant(HH_bin_limits[i]);
        registerParameter(pwa_coefs_reals[i]);
        registerParameter(pwa_coefs_imags[i]);
    }

    resonanceType = 6;
}

__host__ void Spline::recalculateCache() const {
    auto params           = getParameters();
    const unsigned nKnobs = params.size() / 2;
    std::vector<fpcomplex> y(nKnobs);
    unsigned int i = 0;
    for(auto v = params.begin(); v != params.end(); ++v, ++i) {
        unsigned int idx = i / 2;
        fptype value     = parametersList[i];
        if(i % 2 == 0)
            y[idx].real(value);
        else
            y[idx].imag(value);
    }
    //std::vector<fptype> y2_flat = flatten(complex_derivative(constantsList, y));

    //MEMCPY_TO_SYMBOL(cDeriatives, y2_flat.data(), 2 * nKnobs * sizeof(fptype), 0, cudaMemcpyHostToDevice);
}

}

void ResonancePdf::recursiveSetIndices() {
    if(resonanceType == 0) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_RBW");
        GET_FUNCTION_ADDR(ptr_to_RBW);
    } else if(resonanceType == 1) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_GOUSAK");
        GET_FUNCTION_ADDR(ptr_to_GOUSAK);
    } else if(resonanceType == 2) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_LASS");
        GET_FUNCTION_ADDR(ptr_to_LASS);
    } else if(resonanceType == 3) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_NONRES");
        GET_FUNCTION_ADDR(ptr_to_NONRES);
    } else if(resonanceType == 4) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_GAUSSIAN");
        GET_FUNCTION_ADDR(ptr_to_GAUSSIAN);
    }

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

} // namespace GooFit
