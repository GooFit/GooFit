#include <goofit/PDFs/physics/resonances/NonRes_B2KKK.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <complex>
namespace GooFit {

__device__ auto NonRes_B2KKK(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {

    
    fptype alpha            = pc.getParameter(0);
    fptype beta             = pc.getParameter(1);

    fpcomplex ret(0.,0.);
    fpcomplex z1 = exp(-alpha*m12)*fpcomplex(cos(beta*m12),-sin(beta*m12));
    fpcomplex z2 = exp(-alpha*m13)*fpcomplex(cos(beta*m13),-sin(beta*m13));
    
    ret += (z1 + z2);

    pc.incrementIndex(1, 2, 3, 0, 1);
    return ret;
}

__device__ resonance_function_ptr ptr_to_NonRes_B2KKK = NonRes_B2KKK;

namespace Resonances {

NonRes_B2KKK::NonRes_B2KKK(std::string name, Variable ar, Variable ai, Variable alpha, Variable beta)
    : ResonancePdf("NonRes_B2KKK", name, ar, ai) {

    registerParameter(alpha);
    registerParameter(beta);

    registerFunction("ptr_to_NonRes_B2KKK", ptr_to_NonRes_B2KKK);
}

} // namespace Resonances
} // namespace GooFit
