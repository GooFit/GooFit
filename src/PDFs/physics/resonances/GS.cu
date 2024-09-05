#include <goofit/PDFs/physics/resonances/GS.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

namespace GooFit {

template <int I>
__device__ auto gouSak(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);


    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    fptype resmass2 = POW2(resmass);

    fpcomplex result{0.0, 0.0};

    fptype s  = 0.0;
    fptype m = 0.0;
    fptype m1= 0.0;
    fptype m2= 0.0;
    fptype m3= 0.0;


    if(resmass < 1.e-10) {
        GOOFIT_TRACE("Resonance Mass zero!");
        return result;
    }

    if(reswidth < 1.e-10) {
        GOOFIT_TRACE("Resonance Width zero!");
        return result;
    }
    
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
        
        // Implement Gou-Sak
        fptype d_ = d(resmass,q0_);
        fptype f_ = f(m,resmass,reswidth,q_,q0_);

        fptype A = (resmass2-s) + f_;
        fptype B = resmass*gamma;
        fptype C = 1./(A*A + B*B);
        fptype D = 1+(reswidth*d_/resmass);

        fpcomplex ret(A*C,B*C);
        ret *= D;
        ret *= (FR/FR0);
        ret *= (FP/FP0);
        ret *= ZemachSpinFactor;

        result += ret;

        if(I != 0) {
            fptype swpmass = m12;
            m12            = m13;
            m13            = swpmass;
        }
    }

    pc.incrementIndex(1, 2, 2, 0, 1);

    return result;
}

__device__ resonance_function_ptr ptr_to_GOUSAK = gouSak<1>;
__device__ resonance_function_ptr ptr_to_GOUSAK_Sym = gouSak<2>;

namespace Resonances {

GS::GS(std::string name,
       Variable ar,
       Variable ai,
       Variable mass,
       Variable width,
       unsigned int sp,
       unsigned int cyc,
       bool sym)
    : ResonancePdf("GS", name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);

    if(sym)
        registerFunction("ptr_to_GOUSAK_Sym", ptr_to_GOUSAK_Sym);
    else    
        registerFunction("ptr_to_GOUSAK", ptr_to_GOUSAK);
}

} // namespace Resonances
} // namespace GooFit
