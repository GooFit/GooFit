#include <goofit/PDFs/physics/resonances/NonRes3k.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto nonres3k(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    fptype alpha = pc.getParameter(0);

    pc.incrementIndex(1, 1, 0, 0, 1);

    if(sqrt(m12)>2.0 || sqrt(m13)>2.0){
        fptype exp12 = exp(-alpha*m12);
        fptype exp13 = exp(-alpha*m13);
        return fpcomplex(exp12+exp13,0.0);
    }else{
        return fpcomplex(0.0,0.0);
    }
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
