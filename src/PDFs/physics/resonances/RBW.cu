#include <goofit/PDFs/physics/resonances/RBW.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {


template <int I>
__device__ auto plainBW(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);
    bool bachPframe = pc.getConstant(2);
    bool ignoreMom = pc.getConstant(3);
    bool ignoreBW = pc.getConstant(4);
 
    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    fpcomplex result{0.0, 0.0};
    fptype resmass2 = resmass*resmass;

    fptype s  = 0.0;
    fptype m = 0.0;
    fptype m1= 0.0;
    fptype m2= 0.0;
    fptype m3= 0.0;

    // if(resmass < 1.e-10) {
    //     GOOFIT_TRACE("Resonance Mass zero!");
    //     return result;
    // }

    // if(reswidth < 1.e-10) {
    //     GOOFIT_TRACE("Resonance Width zero!");
    //     return result;
    // }

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
            m1 = c_daug3Mass;
            m2 = c_daug1Mass;
            m3 = c_daug2Mass;
        }

         if(PAIR_23 == cyclic_index){
            s    = m23 ;
            m    = sqrt(s);
            m1 = c_daug2Mass;
            m2 = c_daug3Mass;
            m3 = c_daug1Mass;
        }
       

        fptype q0_ = DaugDecayMomResFrame(resmass2, m1, m2);

        fptype p0_ = bachPframe? 
                    BachMomParentFrame(c_motherMass, m3 ,resmass2) :
                    BachMomResFrame(c_motherMass , resmass2, m3);

        fptype q_  = DaugDecayMomResFrame(s, m1, m2);

        fptype p_  =  bachPframe? 
                    BachMomParentFrame(c_motherMass, m3 ,s) : //p*
                    BachMomResFrame(c_motherMass, s , m3); //p

        fptype FR0 = BlattWeisskopfPrime(q0_*c_meson_radius,spin);
        fptype FP0 = BlattWeisskopfPrime(p0_*c_mother_meson_radius,spin);

        fptype FR = BlattWeisskopfPrime(q_*c_meson_radius,spin);
        fptype FP = BlattWeisskopfPrime(p_*c_mother_meson_radius,spin);
        
        fptype gamma = ignoreMom ? reswidth : reswidth*pow(q_/q0_,2.0*spin + 1.0)*(resmass/m)*POW2(FR/FR0);
                  
        fptype A = (-resmass2 + s);
        fptype B = gamma*resmass;
        fptype C = 1.0 / (A*A + B*B);

        fpcomplex ret(A * C, B * C);
        ret *= ignoreBW ? 1.0 : (FR/FR0);
        ret *= ignoreBW ? 1.0 : (FP/FP0);

        fptype cosHel = cFromM(c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m13, m23, m12, cyclic_index);
        fptype legPol = calcLegendrePoly(cosHel,spin);
        fptype ZemachSpinFactor = ignoreMom ? calcZemachSpinFactor(1.0, legPol, spin) :calcZemachSpinFactor(q_*p_, legPol, spin);
        ret *= ZemachSpinFactor;
        //printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f  \t %f  \t %f  \t %f \n",resmass,reswidth,s,m1,m2,m3,c_motherMass,q0_,p0_,q_,p_,FR0,FP0,FR,FP);

        result += ret;

        if(I != 0) {
            fptype swpmass = m13;
            m13            = m23;
            m23            = swpmass;
        }
    }

    
    pc.incrementIndex(1, 2, 5, 0, 1);
  
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
         bool sym,
         bool bachPframe,
         bool ignoreMom,
         bool ignoreBW)
    : ResonancePdf("RBW", name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);
    registerConstant(bachPframe);
    registerConstant(ignoreMom);
    registerConstant(ignoreBW);
 
    if(sym)
        registerFunction("ptr_to_RBW_Sym", ptr_to_RBW_Sym);
    else
        registerFunction("ptr_to_RBW", ptr_to_RBW);
}

} // namespace Resonances

} // namespace GooFit
