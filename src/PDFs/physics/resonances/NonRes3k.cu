#include <goofit/PDFs/physics/resonances/NonRes3k.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto nonres3k(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    fptype alpha = pc.getParameter(0);

    fptype exp12 = exp(-alpha * m12);
    fptype exp13 = exp(-alpha * m13);

    pc.incrementIndex(1, 1, 0, 0, 1);
    return (exp12 + exp13);
}

__device__ resonance_function_ptr ptr_to_NONRES3k = nonres3k;

namespace Resonances {

NonRes3k::NonRes3k(std::string name, Variable ar, Variable ai, Variable alpha)
    : ResonancePdf("NonRes3k", name, ar, ai) {
    registerParameter(alpha);
    registerFunction("ptr_to_NONRES3k", ptr_to_NONRES3k);
}

} // namespace Resonances
} // namespace GooFit
