#include <goofit/PDFs/physics/resonances/RhoOmegaMix.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

template <int I>
__device__ auto rhoomgamix(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);
    bool bachPframe = pc.getConstant(2);
    bool ignoreMom = pc.getConstant(3);
    bool ignoreBW = pc.getConstant(4);
  

    fptype omega_mass  = pc.getParameter(0);
    fptype omega_width = pc.getParameter(1);
    fptype rho_mass    = pc.getParameter(2);
    fptype rho_width   = pc.getParameter(3);

    fpcomplex result{0.0, 0.0};
    fptype omega_mass2 = POW2(omega_mass);
    fptype rho_mass2   = POW2(rho_mass);

    fptype real  = pc.getParameter(4);
    fptype img   = pc.getParameter(5);
    fptype delta = pc.getParameter(6);

    fptype Delta_ = delta * (rho_mass + omega_mass);
    fpcomplex Bterm(real, img);
    Bterm *= Delta_;
    fpcomplex unity(1.0, 0.0);

    if(omega_mass < 1.e-10) {
        GOOFIT_TRACE("Omega Mass zero!");
        return result;
    }

    if(omega_width < 1.e-10) {
        GOOFIT_TRACE("Omega Width zero!");
        return result;
    }

    if(rho_mass < 1.e-10) {
        GOOFIT_TRACE("Rho Mass zero!");
        return result;
    }

    if(rho_width < 1.e-10) {
        GOOFIT_TRACE("Rho Width zero!");
        return result;
    }

#pragma unroll
    for(size_t i = 0; i < I; i++) {
        fptype s    = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
        fptype m          = sqrt(s);
        fptype m1 = PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass;
        fptype m2 = PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass;
        fptype m3
            = (PAIR_12 == cyclic_index ? c_daug3Mass : (PAIR_13 == cyclic_index ? c_daug2Mass : c_daug1Mass));


        // Omega RBW evaluation
        fptype q0_ = DaugDecayMomResFrame(omega_mass2, m1, m2);

        fptype p0_ = bachPframe? 
                        BachMomParentFrame(c_motherMass, m3 ,omega_mass2) :
                        BachMomResFrame(c_motherMass , omega_mass2, m3);
    
        fptype q_  = DaugDecayMomResFrame(s, m1, m2);
    
        fptype p_  =  bachPframe? 
                        BachMomParentFrame(c_motherMass, m3 ,s) : //p*
                        BachMomResFrame(c_motherMass, s , m3); //p
    
        fptype FR0 = BlattWeisskopfPrime(q0_*c_meson_radius,spin);
        fptype FP0 = BlattWeisskopfPrime(p0_*c_mother_meson_radius,spin);
    
        fptype FR = BlattWeisskopfPrime(q_*c_meson_radius,spin);
        fptype FP = BlattWeisskopfPrime(p_*c_mother_meson_radius,spin);
            
        fptype gamma = ignoreMom ? omega_width : omega_width*pow(q_/q0_,2.0*spin + 1.0)*(omega_mass/m)*POW2(FR/FR0);
                      
        fptype A = (omega_mass2 - s);
        fptype B = gamma*omega_mass;
        fptype C = 1.0 / (A*A + B*B);
    
        fpcomplex RBW(A * C, B * C);
        RBW *= ignoreBW ? 1.0 : (FR/FR0);
        RBW *= ignoreBW ? 1.0 : (FP/FP0);

        // end RBW

        // Rho GS evaluation
       
        q0_ = DaugDecayMomResFrame(rho_mass2, m1, m2);

        p0_ = bachPframe? 
                    BachMomParentFrame(c_motherMass, m3 ,rho_mass2) :
                    BachMomResFrame(c_motherMass , rho_mass2, m3);

        q_  = DaugDecayMomResFrame(s, m1, m2);

        p_  =  bachPframe? 
                    BachMomParentFrame(c_motherMass, m3 ,s) : //p*
                    BachMomResFrame(c_motherMass, s , m3); //p

        FR0 = BlattWeisskopfPrime(q0_*c_meson_radius,spin);
        FP0 = BlattWeisskopfPrime(p0_*c_mother_meson_radius,spin);

        FR = BlattWeisskopfPrime(q_*c_meson_radius,spin);
        FP = BlattWeisskopfPrime(p_*c_mother_meson_radius,spin);

        gamma = ignoreMom ? rho_width : rho_width*pow(q_/q0_,2.0*spin + 1.0)*(rho_mass/m)*POW2(FR/FR0);
        
        fptype d_ = d(rho_mass,q0_);
        fptype f_ = f(m,rho_mass,rho_width,q_,q0_);

        A = (rho_mass2-s) + f_;
        B = rho_mass*gamma;
        C = 1./(A*A + B*B);
        fptype D = 1+(rho_width*d_/rho_mass);

        fpcomplex GS(A*C,B*C);
        GS *= D;
        GS *= ignoreBW ? 1.0 : (FR/FR0);
        GS *= ignoreBW ? 1.0 : (FP/FP0);
        // end of Gousak

        // rho-omega mix
        fpcomplex mixingTerm = Bterm * RBW + unity;

        fptype cosHel = cFromM(c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);
        fptype legPol = calcLegendrePoly(cosHel,spin);
        fptype ZemachSpinFactor = calcZemachSpinFactor(q_*p_, legPol, spin);

        result += GS * mixingTerm * ZemachSpinFactor;

        if(I > 1) {
            fptype swpmass = m12;
            m12            = m13;
            m13            = swpmass;
        }
    }
    pc.incrementIndex(1, 7, 5, 0, 1);
    return result;

} // RhoOmegaMix

__device__ resonance_function_ptr ptr_to_RHOOMEGAMIX     = rhoomgamix<1>;
__device__ resonance_function_ptr ptr_to_RHOOMEGAMIX_SYM = rhoomgamix<2>;

namespace Resonances {

RhoOmegaMix::RhoOmegaMix(std::string name,
                         Variable ar,
                         Variable ai,
                         Variable omega_mass,
                         Variable omega_width,
                         Variable rho_mass,
                         Variable rho_width,
                         Variable real,
                         Variable imag,
                         Variable delta,
                         unsigned int sp,
                         unsigned int cyc,
                         bool symDP,
                         bool bachPframe,
                         bool ignoreMom,
                         bool ignoreBW)
    : ResonancePdf("RHOOMEGAMIX", name, ar, ai) {
    registerParameter(omega_mass);
    registerParameter(omega_width);
    registerParameter(rho_mass);
    registerParameter(rho_width);
    registerParameter(real);
    registerParameter(imag);
    registerParameter(delta);

    registerConstant(sp);
    registerConstant(cyc);
    registerConstant(bachPframe);
    registerConstant(ignoreMom);
    registerConstant(ignoreBW);

    if(symDP)
        registerFunction("ptr_to_RHOOMEGAMIX", ptr_to_RHOOMEGAMIX_SYM);
    else
        registerFunction("ptr_to_RHOOMEGAMIX", ptr_to_RHOOMEGAMIX);
}

} // namespace Resonances
} // namespace GooFit
