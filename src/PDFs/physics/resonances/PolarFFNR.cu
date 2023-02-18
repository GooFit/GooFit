#include <goofit/PDFs/physics/resonances/PolarFFNR.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {


template <int I>
__device__ auto func_PolarFFNR(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex {
  
    unsigned int cyclic_index = pc.getConstant(0);
 
    fptype lambda  = pc.getParameter(0);
    
    fptype s  = 0.0;

    fpcomplex ret(0.,0.);

#pragma unroll
    for(size_t i = 0; i < I; i++) {
        
        if(PAIR_12 == cyclic_index){
            s    = m12 ;
        }

        if(PAIR_13 == cyclic_index){
            s    = m13 ;
        }

         if(PAIR_23 == cyclic_index){
            s    = m23 ;
        }
       
        fptype magnitude = 1.0;


        magnitude = 1./(1. + s/(lambda*lambda));

        ret+=fpcomplex(magnitude,0.);

        if(I != 0) {
            fptype swpmass = m13;
            m13            = m23;
            m23            = swpmass;
        }
    }

    
    pc.incrementIndex(1, 1, 1, 0, 1);
  
    return ret;
}

__device__ resonance_function_ptr ptr_to_PolarFFNR     = func_PolarFFNR<1>;
__device__ resonance_function_ptr ptr_to_PolarFFNR_Sym = func_PolarFFNR<2>;

namespace Resonances {

PolarFFNR::PolarFFNR(std::string name,
         Variable ar,
         Variable ai,
         Variable lambda,
         unsigned int cyc,
         bool sym)
    : ResonancePdf("PolarFFNR", name, ar, ai) {
    registerParameter(lambda);
    registerConstant(cyc);

    if(sym)
        registerFunction("ptr_to_PolarFFNR_Sym", ptr_to_PolarFFNR_Sym);
    else
        registerFunction("ptr_to_PolarFFNR", ptr_to_PolarFFNR);
}

} // namespace Resonances

} // namespace GooFit
