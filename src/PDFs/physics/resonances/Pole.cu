#include <goofit/PDFs/physics/resonances/Pole.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {


template <int I>
__device__ auto POLE(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex{
  
    
    fptype real            = pc.getParameter(0);
    fptype img           = pc.getParameter(1);
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
      
        
        fptype reTerm = (real*real - img*img) - s;
        fptype imTerm = 2.0*real*img;

	    fptype scale = 1./(reTerm*reTerm + imTerm*imTerm);
	    fpcomplex ret(reTerm*scale,imTerm*scale);
	
        result += ret;

                                          
        if(I != 0) {
            fptype swpmass = m13;
            m13            = m23;
            m23            = swpmass;
        }
    }

    return result;
}

__device__ resonance_function_ptr ptr_to_POLE      = POLE<1>;
__device__ resonance_function_ptr ptr_to_POLE_SYM  = POLE<2>;

namespace Resonances {

Pole::Pole(std::string name,
         Variable ar,
         Variable ai,
	     Variable real,
         Variable img,
         unsigned int cyc,
         bool symmDP)
    : ResonancePdf("Pole",name, ar, ai) {
    registerParameter(real);
    registerParameter(img);
    registerConstant(cyc);

    if(symmDP) {
        registerFunction("ptr_to_POLE_SYM", ptr_to_POLE_SYM);
    } else {
        registerFunction("ptr_to_POLE", ptr_to_POLE);
    }

}

} // namespace Resonances

} // namespace GooFit
