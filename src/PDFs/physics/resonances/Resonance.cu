#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto dh_dsFun(double s, double daug2Mass, double daug3Mass) -> fptype {
    // Yet another helper function
    const fptype _pi = 3.14159265359;
    double k_s       = twoBodyCMmom(s, daug2Mass, daug3Mass);

    return hFun(s, daug2Mass, daug3Mass) * (1.0 / (8.0 * POW2(k_s)) - 1.0 / (2.0 * s)) + 1.0 / (2.0 * _pi * s);
}

__device__ auto hFun(double s, double daug2Mass, double daug3Mass) -> fptype {
    // Last helper function
    const fptype _pi = 3.14159265359;
    double sm        = daug2Mass + daug3Mass;
    double sqrt_s    = sqrt(s);
    double k_s       = twoBodyCMmom(s, daug2Mass, daug3Mass);

    return ((2 / _pi) * (k_s / sqrt_s) * log((sqrt_s + 2 * k_s) / (sm)));
}

__device__ auto fsFun(double s, double m2, double gam, double daug2Mass, double daug3Mass) -> fptype {
    // Another G-S helper function

    double k_s   = twoBodyCMmom(s, daug2Mass, daug3Mass);
    double k_Am2 = twoBodyCMmom(m2, daug2Mass, daug3Mass);

    double f = gam * m2 / POW3(k_Am2);
    f *= (POW2(k_s) * (hFun(s, daug2Mass, daug3Mass) - hFun(m2, daug2Mass, daug3Mass))
          + (m2 - s) * POW2(k_Am2) * dh_dsFun(m2, daug2Mass, daug3Mass));

    return f;
}

__device__ auto dFun(double s, double daug2Mass, double daug3Mass) -> fptype {
    // Helper function used in Gronau-Sakurai
    const fptype _pi = 3.14159265359;
    double sm        = daug2Mass + daug3Mass;
    double sm24      = sm * sm / 4.0;
    double m         = sqrt(s);
    double k_m2      = twoBodyCMmom(s, daug2Mass, daug3Mass);

    return 3.0 / _pi * sm24 / POW2(k_m2) * log((m + 2 * k_m2) / sm) + m / (2 * _pi * k_m2)
           - sm24 * m / (_pi * POW3(k_m2));
}

__device__ auto twoBodyCMmom(double rMassSq, fptype d1m, fptype d2m) -> fptype {
    // For A -> B + C, calculate momentum of B and C in rest frame of A.
    // PDG 38.16.

    fptype kin1 = 1 - POW2(d1m + d2m) / rMassSq;

    kin1 = kin1 >= 0 ? sqrt(kin1) : 1;

    fptype kin2 = 1 - POW2(d1m - d2m) / rMassSq;
    kin2        = kin2 >= 0 ? sqrt(kin2) : 1;

    return 0.5 * sqrt(rMassSq) * kin1 * kin2;
}

__device__ auto twoBodyCMMothermom(fptype rMassSq, fptype dm, fptype d3m) -> fptype {
    fptype kin1 = 1 - POW2(dm + d3m) / rMassSq;
    if(kin1 >= 0)
        kin1 = sqrt(kin1);
    else
        kin1 = 1;
    fptype kin2 = 1 - POW2(dm - d3m) / rMassSq;
    if(kin2 >= 0)
        kin2 = sqrt(kin2);
    else
        kin2 = 1;

    return 0.5 * rMassSq * kin1 * kin2 / dm;
}

__device__ auto dampingFactorSquare(const fptype &cmmom, const int &spin, const fptype &mRadius) -> fptype {
    fptype square = mRadius * mRadius * cmmom * cmmom;
    fptype dfsq   = 2 * square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    fptype dfsqres = 13 * square / pow(square - 3, 2) + 9 * square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}

__device__ auto dampingFactorSquareNorm(const fptype &cmmom, const int &spin, const fptype &mRadius) -> fptype {
    fptype square = mRadius * mRadius * cmmom * cmmom;
    fptype dfsq   = 1 + square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    fptype dfsqres = dfsq + 8 + 2 * square + square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}

__device__ auto spinFactor(unsigned int spin,
                           fptype motherMass,
                           fptype daug1Mass,
                           fptype daug2Mass,
                           fptype daug3Mass,
                           fptype m12,
                           fptype m13,
                           fptype m23,
                           unsigned int cyclic_index) -> fptype {
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

    fptype massFactor = 1.0 / (_mAB);
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

} // namespace GooFit
