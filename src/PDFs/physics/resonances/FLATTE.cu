#include <goofit/PDFs/physics/resonances/FLATTE.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ void
getAmplitudeCoefficients(fpcomplex a1, fpcomplex a2, fptype &a1sq, fptype &a2sq, fptype &a1a2real, fptype &a1a2imag) {
    // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
    a1sq = thrust::norm(a1);
    a2sq = thrust::norm(a2);
    a1 *= conj(a2);
    a1a2real = a1.real();
    a1a2imag = a1.imag();
}


template <int I>
__device__ auto flatte(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    // indices[1] is unused constant index, for consistency with other function types.
    fptype resmass            = pc.getParameter(0);
    fptype g1                 = pc.getParameter(1);
    fptype g2                 = pc.getParameter(2);
    unsigned int particle     = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);


    const fptype pipmass = 0.13957018;
    const fptype pi0mass = 0.1349766;
    const fptype kpmass  = 0.493677;
    const fptype k0mass  = 0.497614;
    const fptype mEta		= 0.547862; 

    fptype mSumSq0_  = 1.;
    fptype mSumSq1_ = 1.;
    fptype mSumSq2_   = 1.;
    fptype mSumSq3_  = 1.;

    //f0(980)
    if(particle==0){
        mSumSq0_ = 4*pi0mass*pi0mass;
        mSumSq1_ = 4*pipmass*pipmass;
        mSumSq2_ = 4*kpmass*kpmass;
        mSumSq3_ = 4*k0mass*k0mass;
    }
    //a0(980)
    if(particle==1){
        mSumSq0_ = (mEta+pi0mass)*(mEta+pi0mass);
        mSumSq1_ = (mEta+pi0mass)*(mEta+pi0mass);
        mSumSq2_ = (kpmass+kpmass)*(kpmass+kpmass);
        mSumSq3_ = (k0mass+k0mass)*(k0mass+k0mass);
    }

    fpcomplex ret(0., 0.);
    
    fptype rho1(0.0), rho2(0.0);

#pragma unroll    
    for(int i = 0; i < I; i++) {
        fptype s = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
        fptype resmass2 = POW2(resmass);
        fptype dMSq = resmass2 - s;

        if (s > mSumSq0_) {
            rho1 = sqrt(1.0 - mSumSq0_/s)/3.0;
            if (s > mSumSq1_) {
                rho1 += 2.0*sqrt(1.0 - mSumSq1_/s)/3.0;
                if (s > mSumSq2_) {
                    rho2 = 0.5*sqrt(1.0 - mSumSq2_/s);
                    if (s > mSumSq3_) {
                        rho2 += 0.5*sqrt(1.0 - mSumSq3_/s);
                    } else {
                        // Continue analytically below higher channel thresholds
                        // This contributes to the real part of the amplitude denominator
                        dMSq += g2*resmass*0.5*sqrt(mSumSq3_/s - 1.0);
                    }
                } else {
                    // Continue analytically below higher channel thresholds
                    // This contributes to the real part of the amplitude denominator
                    rho2 = 0.0;
                    dMSq += g2*resmass*(0.5*sqrt(mSumSq2_/s - 1.0) + 0.5*sqrt(mSumSq3_/s - 1.0));
                }
            } else {
                // Continue analytically below higher channel thresholds
                // This contributes to the real part of the amplitude denominator
                dMSq += g1*resmass*2.0*sqrt(mSumSq1_/s - 1.0)/3.0;
            }
        }
    
        //the Adler-zero term fA = (m2 − sA)/(m20 − sA) can be used to suppress false 
        //kinematic singularities when m goes below threshold. For f(0)(980), sA = 0.
        
        fptype massFactor = 1.;//resmass;
        
        fptype width1 = g1*rho1*massFactor;
        fptype width2 = g2*rho2*massFactor;
        fptype widthTerm = width1 + width2;
    
        fpcomplex resAmplitude(dMSq, widthTerm);
    
        fptype denomFactor = dMSq*dMSq + widthTerm*widthTerm;
    
        fptype invDenomFactor = 1.0/denomFactor;
    
        resAmplitude *= invDenomFactor;

        ret += resAmplitude;
       
        if(I!=0) {
            fptype swpmass = m12;
            m12            = m13;
            m13            = swpmass;
        }
    }
    pc.incrementIndex(1, 3, 2, 0, 1);
    return ret ;
}

__device__ resonance_function_ptr ptr_to_FLATTE = flatte<1>;
__device__ resonance_function_ptr ptr_to_FLATTE_SYM = flatte<2>;

namespace Resonances {

FLATTE::FLATTE(std::string name,
               Variable ar,
               Variable ai,
               Variable mean,
               Variable g1,
               Variable rg2og1,
               unsigned int particle,
               unsigned int cyc,
               bool symmDP)
    : ResonancePdf("FLATTE", name, ar, ai) {
    registerParameter(mean);
    registerParameter(g1);
    registerParameter(rg2og1);

    registerConstant(particle);
    registerConstant(cyc);
   
    if(symmDP)
        registerFunction("ptr_to_FLATTE_SYM", ptr_to_FLATTE_SYM);
    else
        registerFunction("ptr_to_FLATTE", ptr_to_FLATTE);
}

} // namespace Resonances
} // namespace GooFit
