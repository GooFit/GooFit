#include <goofit/PDFs/physics/resonances/RhoOmegaMix.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

template <int I>
__device__ auto rhoomgamix(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);
    bool norm                 = pc.getConstant(2);

    fptype omega_mass  = pc.getParameter(0);
    fptype omega_width = pc.getParameter(1);
    fptype rho_mass    = pc.getParameter(2);
    fptype rho_width   = pc.getParameter(3);

    fpcomplex result{0.0, 0.0};
    fptype omega_mass2 = POW2(omega_mass);
    fptype rho_mass2   = POW2(rho_mass);

    fptype real  = pc.getParameter(4);
    fptype img   = pc.getParameter(5);
    fptype delta = pc.getParameter(6);

    fptype Delta_ = delta * (rho_mass + omega_mass);
    fpcomplex Bterm(real, img);
    Bterm *= Delta_;
    fpcomplex unity(1.0, 0.0);

#pragma unroll
    for(size_t i = 0; i < I; i++) {
        fptype rMassSq    = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
        fptype m          = sqrt(rMassSq);
        fptype mass_daug1 = PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass;
        fptype mass_daug2 = PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass;
        fptype bachelorMass
            = (PAIR_12 == cyclic_index ? c_daug3Mass : (PAIR_13 == cyclic_index ? c_daug2Mass : c_daug1Mass));

        fptype frFactor       = 1;
        fptype frFactorMother = 1;

        // Calculate momentum of the two daughters in the resonance rest frame
        // Note symmetry under interchange (dm1 <-> dm2)

        fptype measureDaughterMoms = twoBodyCMmom(rMassSq, mass_daug1, mass_daug2);
        fptype nominalDaughterMoms = twoBodyCMmom(omega_mass2, mass_daug1, mass_daug2);

        fptype measureDaughterMomsMother;
        fptype nominalDaughterMomsMother;

        if(norm) {
            // Mother momentum for normalized Blatt-Weisskopf form factors calculated in the resonance rest frame
            measureDaughterMomsMother = twoBodyCMMothermom(rMassSq, c_motherMass, bachelorMass);
            nominalDaughterMomsMother = twoBodyCMMothermom(omega_mass2, c_motherMass, bachelorMass);
        } else {
            // Mother momentum for unnormalized Blatt-Weisskopf form factors calculated in mother rest frame
            measureDaughterMomsMother = twoBodyCMmom(c_motherMass * c_motherMass, sqrt(rMassSq), bachelorMass);
        }

        if(0 != spin) {
            // D0 meson has same spin than resonance
            if(norm) {
                // normalized form factors
                frFactor = dampingFactorSquareNorm(nominalDaughterMoms, spin, c_meson_radius)
                           / dampingFactorSquareNorm(measureDaughterMoms, spin, c_meson_radius);

                frFactorMother = dampingFactorSquareNorm(nominalDaughterMomsMother, spin, c_mother_meson_radius)
                                 / dampingFactorSquareNorm(measureDaughterMomsMother, spin, c_mother_meson_radius);
            }
            // unnormalized form factors
            else {
                frFactor       = dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
                frFactorMother = dampingFactorSquare(measureDaughterMomsMother, spin, c_mother_meson_radius);
            }
        }

        // RBW evaluation
        fptype A = (omega_mass2 - rMassSq);
        fptype B = omega_mass2 * omega_width * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor
                   / sqrt(rMassSq);
        fptype C = 1.0 / (POW2(A) + POW2(B));

        fpcomplex omega(A * C, B * C); // Dropping F_D=1
        omega *= sqrt(frFactor);
        omega *= sqrt(frFactorMother);
        omega *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);
        fptype angular
            = spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

        // Rho GS evaluation
        nominalDaughterMoms = twoBodyCMmom(rho_mass2, mass_daug1, mass_daug2);

        if(norm) {
            // Mother momentum for normalized Blatt-Weisskopf form factors calculated in the resonance rest frame
            measureDaughterMomsMother = twoBodyCMMothermom(rMassSq, c_motherMass, bachelorMass);
            nominalDaughterMomsMother = twoBodyCMMothermom(rho_mass2, c_motherMass, bachelorMass);
        } else {
            // Mother momentum for unnormalized Blatt-Weisskopf form factors calculated in mother rest frame
            measureDaughterMomsMother = twoBodyCMmom(c_motherMass * c_motherMass, sqrt(rMassSq), bachelorMass);
        }
        if(0 != spin) {
            if(norm) {
                frFactor = dampingFactorSquareNorm(nominalDaughterMoms, spin, c_meson_radius);
                frFactor /= dampingFactorSquareNorm(measureDaughterMoms, spin, c_meson_radius);

                frFactorMother = dampingFactorSquareNorm(nominalDaughterMomsMother, spin, c_mother_meson_radius);
                frFactorMother /= dampingFactorSquareNorm(measureDaughterMomsMother, spin, c_mother_meson_radius);
            }
            // unnormalized form factors
            else {
                frFactor       = dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
                frFactorMother = dampingFactorSquare(measureDaughterMomsMother, spin, c_mother_meson_radius);
            }
        }
        // Implement Gou-Sak:

        fptype D = (1.0 + dFun(rho_mass2, c_daug2Mass, c_daug3Mass) * rho_width / sqrt(rho_mass2));
        fptype E = rho_mass2 - rMassSq + fsFun(rMassSq, rho_mass2, rho_width, c_daug2Mass, c_daug3Mass);
        fptype F
            = sqrt(rho_mass2) * rho_width * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor;

        D /= (E * E + F * F);
        fpcomplex rho(D * E, D * F); // Dropping F_D=1
        rho *= sqrt(frFactor);
        rho *= sqrt(frFactorMother);
        rho *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);
        // end of Gousak

        // rho-omega mix
        fpcomplex mixingTerm = Bterm * omega + unity;
        result += rho * mixingTerm;
    }
    pc.incrementIndex(1, 7, 3, 0, 1);
    return result;

} // RhoOmegaMix

__device__ resonance_function_ptr ptr_to_RHOOMEGAMIX = rhoomgamix<1>;

namespace Resonances {

RhoOmegaMix::RhoOmegaMix(std::string name,
                         Variable ar,
                         Variable ai,
                         Variable omega_mass,
                         Variable omega_width,
                         Variable rho_mass,
                         Variable rho_width,
                         Variable real,
                         Variable imag,
                         Variable delta,
                         unsigned int sp,
                         unsigned int cyc,
                         bool norm,
                         bool sym)
    : ResonancePdf("RHOOMEGAMIX", name, ar, ai) {
    registerParameter(omega_mass);
    registerParameter(omega_width);
    registerParameter(rho_mass);
    registerParameter(rho_width);
    registerParameter(real);
    registerParameter(imag);
    registerParameter(delta);

    registerConstant(sp);
    registerConstant(cyc);

    registerConstant(norm);

    if(sym)
        registerFunction("ptr_to_RHOOMEGAMIX", ptr_to_RHOOMEGAMIX);
    else
        registerFunction("ptr_to_RHOOMEGAMIX", ptr_to_RHOOMEGAMIX);
}

} // namespace Resonances
} // namespace GooFit
