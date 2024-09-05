#include <goofit/PDFs/physics/resonances/NonRes3k.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>
namespace GooFit {

__device__ auto nonres3k(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex {
    fptype alpha = pc.getParameter(0);
    fptype beta = pc.getParameter(1);

    pc.incrementIndex(1, 2, 0, 0, 1);

    // if(m23 > 4.0 && m13> 4.0 ) {
        fptype exp13= exp(-alpha * m13);
        fptype exp23 = exp(-alpha * m23);
        fpcomplex amp_m13 = exp13*fpcomplex(cos(-beta*m13),sin(-beta*m13));
        fpcomplex amp_m23 = exp23*fpcomplex(cos(-beta*m23),sin(-beta*m23));
        return amp_m13+amp_m23;
    // } else {
    //    return fpcomplex(0.00001, 0.0);
    // }
}

__device__ resonance_function_ptr ptr_to_NONRES3k = nonres3k;

namespace Resonances {

NonRes3k::NonRes3k(std::string name, Variable ar, Variable ai, Variable alpha, Variable beta)
    : ResonancePdf("NonRes3k", name, ar, ai) {
    registerParameter(alpha);
    registerParameter(beta);
    registerFunction("ptr_to_NONRES3k", ptr_to_NONRES3k);
}

} // namespace Resonances
} // namespace GooFit
