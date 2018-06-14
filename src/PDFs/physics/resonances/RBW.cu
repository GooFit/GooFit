#include <goofit/PDFs/physics/resonances/RBW.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

template <int I>
__device__ fpcomplex plainBW(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);

    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    fpcomplex result{0.0, 0.0};
    fptype resmass2 = POW2(resmass);

#pragma unroll
    for(size_t i = 0; i < I; i++) {
        fptype rMassSq    = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
        fptype mass_daug1 = PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass;
        fptype mass_daug2 = PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass;

        fptype frFactor = 1;

        // Calculate momentum of the two daughters in the resonance rest frame
        // Note symmetry under interchange (dm1 <-> dm2)

        fptype measureDaughterMoms = twoBodyCMmom(rMassSq, mass_daug1, mass_daug2);
        fptype nominalDaughterMoms = twoBodyCMmom(resmass2, mass_daug1, mass_daug2);

        if(0 != spin) {
            frFactor = dampingFactorSquare(nominalDaughterMoms, spin, c_meson_radius)
                       / dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
        }

        // RBW evaluation
        fptype A = (resmass2 - rMassSq);
        fptype B = resmass2 * reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor
                   / sqrt(rMassSq);
        fptype C = 1.0 / (POW2(A) + POW2(B));

        fpcomplex ret(A * C, B * C); // Dropping F_D=1

        ret *= sqrt(frFactor);
        ret *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

        result += ret;

        if(I > 1) {
            cyclic_index = cyclic_index + 1 % 3;
        }
    }
    pc.incrementIndex(1, 2, 2, 0, 1);
    return result;
}

__device__ resonance_function_ptr ptr_to_RBW     = plainBW<1>;
__device__ resonance_function_ptr ptr_to_RBW_Sym = plainBW<2>;

namespace Resonances {

RBW::RBW(std::string name,
         Variable ar,
         Variable ai,
         Variable mass,
         Variable width,
         unsigned int sp,
         unsigned int cyc,
         bool sym)
    : ResonancePdf(name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);

    if(sym)
        registerFunction("ptr_to_RBW_Sym", ptr_to_RBW_Sym);
    else
        registerFunction("ptr_to_RBW", ptr_to_RBW);
}

} // namespace Resonances

} // namespace GooFit
