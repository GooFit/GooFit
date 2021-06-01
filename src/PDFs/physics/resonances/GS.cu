#include <goofit/PDFs/physics/resonances/GS.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ fptype hFun(double s, double daug2Mass, double daug3Mass) {
    // Last helper function
    const fptype _pi = 3.14159265359;
    double sm        = daug2Mass + daug3Mass;
    double sqrt_s    = sqrt(s);
    double k_s       = twoBodyCMmom(s, daug2Mass, daug3Mass);

    return ((2 / _pi) * (k_s / sqrt_s) * log((sqrt_s + 2 * k_s) / (sm)));
}

__device__ fptype dh_dsFun(double s, double daug2Mass, double daug3Mass) {
    // Yet another helper function
    const fptype _pi = 3.14159265359;
    double k_s       = twoBodyCMmom(s, daug2Mass, daug3Mass);

    return hFun(s, daug2Mass, daug3Mass) * (1.0 / (8.0 * POW2(k_s)) - 1.0 / (2.0 * s)) + 1.0 / (2.0 * _pi * s);
}

__device__ fptype dFun(double s, double daug2Mass, double daug3Mass) {
    // Helper function used in Gronau-Sakurai
    const fptype _pi = 3.14159265359;
    double sm        = daug2Mass + daug3Mass;
    double sm24      = sm * sm / 4.0;
    double m         = sqrt(s);
    double k_m2      = twoBodyCMmom(s, daug2Mass, daug3Mass);

    return 3.0 / _pi * sm24 / POW2(k_m2) * log((m + 2 * k_m2) / sm) + m / (2 * _pi * k_m2)
           - sm24 * m / (_pi * POW3(k_m2));
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
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);
    bool norm                 = pc.getConstant(2);

    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype bachelorMass
        = (PAIR_12 == cyclic_index ? c_daug3Mass : (PAIR_13 == cyclic_index ? c_daug2Mass : c_daug1Mass));
    fptype frFactor  = 1;
    fptype frFactorD = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).
    fptype measureDaughterMoms = twoBodyCMmom(rMassSq,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(resmass,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass));

    fptype measureDaughterMomsMother;
    fptype nominalDaughterMomsMother;

    if(norm) {
        // Mother momentum for normalized Blatt-Weisskopf form factors calculated in the resonance rest frame
        measureDaughterMomsMother = twoBodyCMMothermom(rMassSq, c_motherMass, bachelorMass);
        nominalDaughterMomsMother = twoBodyCMMothermom(resmass, c_motherMass, bachelorMass);
    } else {
        // Mother momentum for unnormalized Blatt-Weisskopf form factors calculated in mother rest frame
        measureDaughterMomsMother = twoBodyCMmom(c_motherMass * c_motherMass, sqrt(rMassSq), bachelorMass);
    }

    if(0 != spin) {
        if(norm) {
            frFactor = dampingFactorSquareNorm(nominalDaughterMoms, spin, c_meson_radius);
            frFactor /= dampingFactorSquareNorm(measureDaughterMoms, spin, c_meson_radius);

            frFactorD = dampingFactorSquareNorm(nominalDaughterMomsMother, spin, c_mother_meson_radius);
            frFactorD /= dampingFactorSquareNorm(measureDaughterMomsMother, spin, c_mother_meson_radius);
        }
        // unnormalized form factors
        else {
            frFactor  = dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
            frFactorD = dampingFactorSquare(measureDaughterMomsMother, spin, c_mother_meson_radius);
        }
    }

    // Implement Gou-Sak:

    fptype D = (1.0 + dFun(resmass, c_daug2Mass, c_daug3Mass) * reswidth / sqrt(resmass));
    fptype E = resmass - rMassSq + fsFun(rMassSq, resmass, reswidth, c_daug2Mass, c_daug3Mass);
    fptype F = sqrt(resmass) * reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor;

    D /= (E * E + F * F);
    fpcomplex retur(D * E, D * F); // Dropping F_D=1
    retur *= sqrt(frFactor);
    retur *= sqrt(frFactorD);
    retur *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

    pc.incrementIndex(1, 2, 3, 0, 1);

    return retur;
}

__device__ resonance_function_ptr ptr_to_GOUSAK = gouSak;

namespace Resonances {

GS::GS(std::string name,
       Variable ar,
       Variable ai,
       Variable mass,
       Variable width,
       unsigned int sp,
       unsigned int cyc,
       bool norm,
       bool sym)
    : ResonancePdf("GS", name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);
    registerConstant(norm);

    registerFunction("ptr_to_GOUSAK", ptr_to_GOUSAK);
}

} // namespace Resonances
} // namespace GooFit
