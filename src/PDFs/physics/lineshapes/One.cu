#include <goofit/PDFs/physics/lineshapes/One.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto LS_ONE(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    pc.incrementIndex(1, 0, 1, 0, 1);
    return {1., 0.};
}

__device__ resonance_function_ptr ptr_to_LS_ONE = LS_ONE;

Lineshapes::One::One(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("One", name, L, Mpair, FormFac, radius) {
    // TODO: Simplify signature!

    registerFunction("ptr_to_LS_ONE", ptr_to_LS_ONE);

    initialize();
}

bool Lineshapes::One::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

} // namespace GooFit
