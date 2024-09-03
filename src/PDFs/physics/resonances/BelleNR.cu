#include <goofit/PDFs/physics/resonances/Pole.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {


template <int I>
__device__ auto BELLENR(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex{
  
    
    fptype alpha            = pc.getParameter(0);
    unsigned int cyclic_index = pc.getConstant(0);
    
    fpcomplex result(0., 0.);
    fptype s  = 0.0;

#pragma unroll
    for(int i = 0; i < I; i++) {
        if(PAIR_12 == cyclic_index){
            s    = m12 ;
        }

        if(PAIR_13 == cyclic_index){
            s    = m13 ;
        }
         if(PAIR_23 == cyclic_index){
            s    = m23 ;
        }
	
        result += fpcomplex(exp(-alpha*s),0.0);

                                          
        if(I != 0) {
            fptype swpmass = m13;
            m13            = m23;
            m23            = swpmass;
        }
    }

    return result;
}

__device__ resonance_function_ptr ptr_to_BELLENR      = BELLENR<1>;
__device__ resonance_function_ptr ptr_to_BELLENR_SYM  = BELLENR<2>;

namespace Resonances {

BelleNR::BelleNR(std::string name,
         Variable ar,
         Variable ai,
	     Variable alpha,
         unsigned int cyc,
         bool symmDP)
    : ResonancePdf("BelleNR",name, ar, ai) {
    registerParameter(alpha);
    registerConstant(cyc);

    if(symmDP) {
        registerFunction("ptr_to_BELLENR_SYM", ptr_to_BELLENR_SYM);
    } else {
        registerFunction("ptr_to_BELLENR", ptr_to_BELLENR);
    }

}

} // namespace Resonances

} // namespace GooFit
