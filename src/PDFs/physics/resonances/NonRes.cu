#include <goofit/PDFs/physics/resonances/NonRes.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto nonres(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    pc.incrementIndex(1, 0, 0, 0, 1);
    return {1.0, 0.0};
}

__device__ resonance_function_ptr ptr_to_NONRES = nonres;

namespace Resonances {

NonRes::NonRes(std::string name, Variable ar, Variable ai)
    : ResonancePdf("NonRes", name, ar, ai) {
    registerFunction("ptr_to_NONRES", ptr_to_NONRES);
}

} // namespace Resonances
} // namespace GooFit
