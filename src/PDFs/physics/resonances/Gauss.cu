#include <goofit/PDFs/physics/resonances/Gauss.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto gaussian(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    // indices[1] is unused constant index, for consistency with other function types.
    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    unsigned int cyclic_index = pc.getConstant(0);

    // Notice sqrt - this function uses mass, not mass-squared like the other resonance types.
    fptype massToUse = sqrt(PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    massToUse -= resmass;
    massToUse /= reswidth;
    massToUse *= massToUse;
    fptype ret = exp(-0.5 * massToUse);

    pc.incrementIndex(1, 2, 1, 0, 1);
    // Ignore factor 1/sqrt(2pi).
    ret /= reswidth;

    return {ret, 0.0};
}

__device__ resonance_function_ptr ptr_to_GAUSSIAN = gaussian;

namespace Resonances {

// Constructor for regular BW,Gounaris-Sakurai,LASS
Gauss::Gauss(std::string name, Variable ar, Variable ai, Variable mass, Variable width, unsigned int cyc)
    : ResonancePdf("Gauss", name, ar, ai) {
    // Making room for index of decay-related constants. Assumption:
    // These are mother mass and three daughter masses in that order.
    // They will be registered by the object that uses this resonance,
    // which will tell this object where to find them by calling setConstantIndex.
    registerParameter(mass);
    registerParameter(width);

    registerConstant(cyc);

    registerFunction("ptr_to_GAUSSIAN", ptr_to_GAUSSIAN);
}

} // namespace Resonances
} // namespace GooFit
