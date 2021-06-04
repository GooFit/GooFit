#include <goofit/PDFs/physics/resonances/LASS.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto lass(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);
    bool norm                 = pc.getConstant(2);

    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    fptype _a    = pc.getParameter(2);
    fptype _r    = pc.getParameter(3);
    fptype _R    = pc.getParameter(4);
    fptype _phiR = pc.getParameter(5);
    fptype _B    = pc.getParameter(6);
    fptype _phiB = pc.getParameter(7);

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
                                              (PAIR_23 == cyclic_index ? c_daug3Mass : c_daug2Mass));
    fptype nominalDaughterMoms = twoBodyCMmom(resmass,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_23 == cyclic_index ? c_daug3Mass : c_daug2Mass));

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
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, c_meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);

        frFactorD = dampingFactorSquare(nominalDaughterMomsMother, spin, c_mother_meson_radius);
        frFactorD /= dampingFactorSquare(measureDaughterMomsMother, spin, c_mother_meson_radius);
        if(norm) {
            // normalized form factors
            frFactor = dampingFactorSquareNorm(nominalDaughterMoms, spin, c_meson_radius)
                       / dampingFactorSquareNorm(measureDaughterMoms, spin, c_meson_radius);

            frFactorD = dampingFactorSquareNorm(nominalDaughterMomsMother, spin, c_mother_meson_radius)
                        / dampingFactorSquareNorm(measureDaughterMomsMother, spin, c_mother_meson_radius);
        }
        // unnormalized form factors
        else {
            frFactor  = dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
            frFactorD = dampingFactorSquare(measureDaughterMomsMother, spin, c_mother_meson_radius);
        }
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

    /*
       fptype _a    = 0.113; //0.22357;
       fptype _r    = -33.8; //-15.042;
       fptype _R    = 1; // ? Fixed
       fptype _phiR = -1.91463; //1.10644;
       fptype _B    = 0.96; //0.614463;
       fptype _phiB = 0.00174533; //-0.0981907;
     */

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
    resT *= sqrt(frFactorD);
    resT *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

    pc.incrementIndex(1, 8, 3, 0, 1);

    return resT;
}

__device__ resonance_function_ptr ptr_to_LASS = lass;

namespace Resonances {

LASS::LASS(std::string name,
           Variable ar,
           Variable ai,
           Variable mass,
           Variable width,
           Variable _a,
           Variable _r,
           Variable _R,
           Variable _phiR,
           Variable _B,
           Variable _phiB,
           unsigned int sp,
           unsigned int cyc,
           bool norm)
    : ResonancePdf("LASS", name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);
    registerParameter(_a);
    registerParameter(_r);
    registerParameter(_R);
    registerParameter(_phiR);
    registerParameter(_B);
    registerParameter(_phiB);

    registerConstant(sp);
    registerConstant(cyc);
    registerConstant(norm);

    registerFunction("ptr_to_LASS", ptr_to_LASS);
}

} // namespace Resonances
} // namespace GooFit
