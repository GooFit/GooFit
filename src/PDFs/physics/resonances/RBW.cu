#include <goofit/PDFs/physics/resonances/RBW.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>
namespace GooFit {

// We use here the covariant formalism for spin: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.51.2247

template <int I>
__device__ auto plainBW(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex {
    fptype spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);
 
    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    fpcomplex result{0.0, 0.0};
   
    fptype s  = 0.0;
    fptype m = 0.0;
    fptype m1= 0.0;
    fptype m2= 0.0;
    fptype m3= 0.0;

    fptype resmass2 = resmass*resmass;
    
#pragma unroll
    for(size_t i = 0; i < I; i++) {
        
        if(PAIR_12 == cyclic_index){
            s    = m12 ;
            m    = sqrt(s);
            m1 = c_daug1Mass;
            m2 = c_daug2Mass;
            m3 = c_daug3Mass;
        }

        if(PAIR_13 == cyclic_index){
            s    = m13 ;
            m    = sqrt(s);
            m1 = c_daug1Mass;
            m2 = c_daug3Mass;
            m3 = c_daug2Mass;
        }
         if(PAIR_23 == cyclic_index){
            s    = m23 ;
            m    = sqrt(s);
            m1 = c_daug2Mass;
            m2 = c_daug3Mass;
            m3 = c_daug1Mass;
        }

        fptype erm = getERM(c_motherMass, m3, m );
        fptype erm0 = getERM(c_motherMass, m3, resmass );
        fptype covFactor0 = calcCovFactor( erm0 , spin);


        fptype q0_ = DaugDecayMomResFrame(resmass2, m1, m2);
        fptype p0_ =  BachMomParentFrame(c_motherMass, m3 ,resmass2) ; //p0*

        if(spin>=2.)
            covFactor0 = pow(covFactor0, 1.0/spin );  

        fptype FR0 = BlattWeisskopfPrime(q0_*c_meson_radius,spin);
        fptype FP0 = BlattWeisskopfPrime(p0_*covFactor0*c_mother_meson_radius,spin);

        fptype q_ = calc_q(s, m2, m1);
        fptype p_  =  BachMomParentFrame(c_motherMass, m3 ,s) ; //p*
        fptype covFactor = calcCovFactor( erm , spin);

        fptype cosHel = cFromM(c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m13, m23, m12, cyclic_index);
        fptype legPol = calcLegendrePoly(cosHel,spin);
        fptype ZemachSpinFactor = calcZemachSpinFactor(q_*p_, legPol, spin);
        
        if(spin>=2.)
            covFactor = pow(covFactor, 1.0/spin );  

        ZemachSpinFactor *= covFactor;
        
        fptype FR = BlattWeisskopfPrime(q_*c_meson_radius,spin);
        fptype FP = BlattWeisskopfPrime(p_*covFactor*c_mother_meson_radius,spin);
        fptype gamma = reswidth*pow(q_/q0_,2.0*spin + 1.0)*(resmass/m)*POW2(FR/FR0);
                  
        fptype A = (resmass2 - s);
        fptype B = gamma*resmass;
        fptype C = 1.0 / (A*A + B*B);

        fpcomplex ret(A * C, B * C);
        ret *=  (FR/FR0);
        ret *=  (FP/FP0);
        ret *= ZemachSpinFactor;
    
        result += ret;

        if(I != 0) {
            fptype swpmass = m13;
            m13            = m23;
            m23            = swpmass;
        }
    }

    pc.incrementIndex(1, 2, 2, 0, 1);
  
    return result;
}

__device__ resonance_function_ptr ptr_to_RBW     = plainBW<1>;
__device__ resonance_function_ptr ptr_to_RBW_Sym = plainBW<2>;

namespace Resonances {

RBW::RBW(std::string name,
         Variable ar,
         Variable ai,
         Variable mass,
         Variable width,
         unsigned int sp,
         unsigned int cyc,
         bool sym)
    : ResonancePdf("RBW", name, ar, ai) {

    registerParameter(mass);
    registerParameter(width);
    registerConstant(sp);
    registerConstant(cyc);
 
    if(sym)
        registerFunction("ptr_to_RBW_Sym", ptr_to_RBW_Sym);
    else
        registerFunction("ptr_to_RBW", ptr_to_RBW);
}

} // namespace Resonances

} // namespace GooFit
