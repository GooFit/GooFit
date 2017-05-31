#include "goofit/PDFs/physics/ResonancePdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

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

__device__ thrust::complex<fptype> plainBW(fptype m12, fptype m13, fptype m23, unsigned int *indices) {
    fptype motherMass   = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 0]);
    fptype daug1Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 1]);
    fptype daug2Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 2]);
    fptype daug3Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 3]);
    fptype meson_radius = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 4]);

    fptype resmass            = RO_CACHE(cudaArray[RO_CACHE(indices[2])]);
    fptype reswidth           = RO_CACHE(cudaArray[RO_CACHE(indices[3])]);
    unsigned int spin         = RO_CACHE(indices[4]);
    unsigned int cyclic_index = RO_CACHE(indices[5]);

    fptype rMassSq  = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype frFactor = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).
    fptype measureDaughterMoms = twoBodyCMmom(
        rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(
        resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

    if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius);
    }

    // RBW evaluation
    fptype A = (resmass - rMassSq);
    fptype B = resmass * reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor
               / sqrt(rMassSq);
    fptype C = 1.0 / (A * A + B * B);
    thrust::complex<fptype> ret(A * C, B * C); // Dropping F_D=1

    ret *= sqrt(frFactor);
    fptype spinF = spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);
    ret *= spinF;
    // printf("%f, %f, %f, %f\n",ret.real, ret.imag, m12, m13);
    return ret;
}

__device__ thrust::complex<fptype> gaussian(fptype m12, fptype m13, fptype m23, unsigned int *indices) {
    // indices[1] is unused constant index, for consistency with other function types.
    fptype resmass            = RO_CACHE(cudaArray[RO_CACHE(indices[2])]);
    fptype reswidth           = RO_CACHE(cudaArray[RO_CACHE(indices[3])]);
    unsigned int cyclic_index = indices[4];

    // Notice sqrt - this function uses mass, not mass-squared like the other resonance types.
    fptype massToUse = sqrt(PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    massToUse -= resmass;
    massToUse /= reswidth;
    massToUse *= massToUse;
    fptype ret = exp(-0.5 * massToUse);

    // Ignore factor 1/sqrt(2pi).
    ret /= reswidth;

    return thrust::complex<fptype>(ret, 0);
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

__device__ thrust::complex<fptype> gouSak(fptype m12, fptype m13, fptype m23, unsigned int *indices) {
    fptype motherMass   = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 0]);
    fptype daug1Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 1]);
    fptype daug2Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 2]);
    fptype daug3Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 3]);
    fptype meson_radius = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 4]);

    fptype resmass            = RO_CACHE(cudaArray[RO_CACHE(indices[2])]);
    fptype reswidth           = RO_CACHE(cudaArray[RO_CACHE(indices[3])]);
    unsigned int spin         = RO_CACHE(indices[4]);
    unsigned int cyclic_index = RO_CACHE(indices[5]);

    fptype rMassSq  = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype frFactor = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).
    fptype measureDaughterMoms = twoBodyCMmom(
        rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(
        resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

    if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius);
    }

    // Implement Gou-Sak:

    fptype D = (1.0 + dFun(resmass, daug2Mass, daug3Mass) * reswidth / sqrt(resmass));
    fptype E = resmass - rMassSq + fsFun(rMassSq, resmass, reswidth, daug2Mass, daug3Mass);
    fptype F = sqrt(resmass) * reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor;

    D /= (E * E + F * F);
    thrust::complex<fptype> retur(D * E, D * F); // Dropping F_D=1
    retur *= sqrt(frFactor);
    retur *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);

    return retur;
}

__device__ thrust::complex<fptype> lass(fptype m12, fptype m13, fptype m23, unsigned int *indices) {
    fptype motherMass   = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 0]);
    fptype daug1Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 1]);
    fptype daug2Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 2]);
    fptype daug3Mass    = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 3]);
    fptype meson_radius = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 4]);

    fptype resmass            = RO_CACHE(cudaArray[RO_CACHE(indices[2])]);
    fptype reswidth           = RO_CACHE(cudaArray[RO_CACHE(indices[3])]);
    unsigned int spin         = RO_CACHE(indices[4]);
    unsigned int cyclic_index = RO_CACHE(indices[5]);

    fptype rMassSq  = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype frFactor = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).

    fptype measureDaughterMoms = twoBodyCMmom(
        rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_23 == cyclic_index ? daug3Mass : daug2Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(
        resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_23 == cyclic_index ? daug3Mass : daug2Mass));

    if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius);
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
    thrust::complex<fptype> expi2deltaB
        = thrust::complex<fptype>(qcot_deltaB, q) / thrust::complex<fptype>(qcot_deltaB, -q);
    thrust::complex<fptype> resT = thrust::complex<fptype>(cos(_phiR + 2 * _phiB), sin(_phiR + 2 * _phiB)) * _R;

    thrust::complex<fptype> prop
        = thrust::complex<fptype>(1, 0) / thrust::complex<fptype>(resmass - rMassSq, sqrt(resmass) * g);
    // resT *= prop*m0*_g0*m0/twoBodyCMmom(m0*m0, _trackinfo[i])*expi2deltaB;
    resT *= prop * (resmass * reswidth / nominalDaughterMoms) * expi2deltaB;

    // calculate bkg part
    resT += thrust::complex<fptype>(cos(_phiB), sin(_phiB)) * _B * (cos(_phiB) + cot_deltaB * sin(_phiB))
            * sqrt(rMassSq) / thrust::complex<fptype>(qcot_deltaB, -q);

    resT *= sqrt(frFactor);
    resT *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);

    return resT;
}

__device__ thrust::complex<fptype> nonres(fptype m12, fptype m13, fptype m23, unsigned int *indices) {
    return thrust::complex<fptype>(1, 0);
}

__device__ void getAmplitudeCoefficients(thrust::complex<fptype> a1,
                                         thrust::complex<fptype> a2,
                                         fptype &a1sq,
                                         fptype &a2sq,
                                         fptype &a1a2real,
                                         fptype &a1a2imag) {
    // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
    a1sq = thrust::norm(a1);
    a2sq = thrust::norm(a2);
    a1 *= conj(a2);
    a1a2real = a1.real();
    a1a2imag = a1.imag();
}

__device__ resonance_function_ptr ptr_to_RBW      = plainBW;
__device__ resonance_function_ptr ptr_to_GOUSAK   = gouSak;
__device__ resonance_function_ptr ptr_to_GAUSSIAN = gaussian;
__device__ resonance_function_ptr ptr_to_NONRES   = nonres;
__device__ resonance_function_ptr ptr_to_LASS     = lass;

ResonancePdf::ResonancePdf(
    std::string name, Variable *ar, Variable *ai, Variable *mass, Variable *width, unsigned int sp, unsigned int cyc)
    : GooPdf(nullptr, name)
    , amp_real(ar)
    , amp_imag(ai) {
    std::vector<unsigned int> pindices;
    pindices.push_back(0);
    // Making room for index of decay-related constants. Assumption:
    // These are mother mass and three daughter masses in that order.
    // They will be registered by the object that uses this resonance,
    // which will tell this object where to find them by calling setConstantIndex.

    pindices.push_back(registerParameter(mass));
    pindices.push_back(registerParameter(width));
    pindices.push_back(sp);
    pindices.push_back(cyc);

    GET_FUNCTION_ADDR(ptr_to_RBW);
    initialize(pindices);
}

ResonancePdf::ResonancePdf(
    std::string name, Variable *ar, Variable *ai, unsigned int sp, Variable *mass, Variable *width, unsigned int cyc)
    : GooPdf(nullptr, name)
    , amp_real(ar)
    , amp_imag(ai) {
    // Same as BW except for function pointed to.
    std::vector<unsigned int> pindices;
    pindices.push_back(0);
    pindices.push_back(registerParameter(mass));
    pindices.push_back(registerParameter(width));
    pindices.push_back(sp);
    pindices.push_back(cyc);

    GET_FUNCTION_ADDR(ptr_to_GOUSAK);
    initialize(pindices);
}

ResonancePdf::ResonancePdf(
    std::string name, Variable *ar, Variable *ai, Variable *mass, unsigned int sp, Variable *width, unsigned int cyc)
    : GooPdf(nullptr, name)
    , amp_real(ar)
    , amp_imag(ai) {
    // Same as BW except for function pointed to.
    std::vector<unsigned int> pindices;
    pindices.push_back(0);
    pindices.push_back(registerParameter(mass));
    pindices.push_back(registerParameter(width));
    pindices.push_back(sp);
    pindices.push_back(cyc);

    GET_FUNCTION_ADDR(ptr_to_LASS);
    initialize(pindices);
}

ResonancePdf::ResonancePdf(std::string name, Variable *ar, Variable *ai)
    : GooPdf(nullptr, name)
    , amp_real(ar)
    , amp_imag(ai) {
    std::vector<unsigned int> pindices;
    pindices.push_back(0);
    // Dummy index for constants - won't use it, but calling
    // functions can't know that and will call setConstantIndex anyway.
    GET_FUNCTION_ADDR(ptr_to_NONRES);
    initialize(pindices);
}

ResonancePdf::ResonancePdf(
    std::string name, Variable *ar, Variable *ai, Variable *mean, Variable *sigma, unsigned int cyc)
    : GooPdf(nullptr, name)
    , amp_real(ar)
    , amp_imag(ai) {
    std::vector<unsigned int> pindices;
    pindices.push_back(0);
    // Dummy index for constants - won't use it, but calling
    // functions can't know that and will call setConstantIndex anyway.
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigma));
    pindices.push_back(cyc);

    GET_FUNCTION_ADDR(ptr_to_GAUSSIAN);

    initialize(pindices);
}

} // namespace GooFit
