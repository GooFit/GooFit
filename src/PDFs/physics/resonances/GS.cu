#include <goofit/PDFs/physics/resonances/GS.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

template <int I>
__device__ auto gouSak(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);
    bool norm                 = pc.getConstant(2);

    fpcomplex ret(0.,0.);

    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    
    if(resmass < 1.e-10){
        GOOFIT_TRACE("Resonance Mass zero!");
        return ret;
    }

    if(reswidth < 1.e-10){
        GOOFIT_TRACE("Resonance Width zero!");
        return ret;
    }

    #pragma unroll
    for(size_t i = 0; i < I; i++) {

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

        ret += retur;
        
        if(I > 1) {
             fptype swpmass = m12;
            m12            = m13;
            m13            = swpmass;
        }
    }

    pc.incrementIndex(1, 2, 3, 0, 1);

    return ret;
}

__device__ resonance_function_ptr ptr_to_GOUSAK = gouSak<1>;
__device__ resonance_function_ptr ptr_to_GOUSAK_SYM = gouSak<2>;

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

    if(sym)
        registerFunction("ptr_to_GOUSAK", ptr_to_GOUSAK_SYM);
    else
        registerFunction("ptr_to_GOUSAK", ptr_to_GOUSAK);
}

} // namespace Resonances
} // namespace GooFit
