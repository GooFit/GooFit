#include <goofit/PDFs/physics/resonances/ScatteringAmp.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {
__device__ fptype Deriatives[2 * 100];
__device__ fptype phi_s[2*100];
__device__ fptype phi_real[2*100];
__device__ fptype phi_imag[2*100];

///////////Amp PIPI->PIPI equation 7b (PHYSICAL REVIEW LETTERS 131, 051802 (2023))

/// Flatten a complex array into a standard one (1r, 1i, 2r, 2i, ...)
template <typename T>
auto make_flatten(const std::vector<thrust::complex<T>> &input) -> std::vector<T> {
    std::vector<T> output;
    for(auto val : input) {
        output.push_back(val.real());
        output.push_back(val.imag());
    }
    return output;
}

auto calc_complex_derivative(const std::vector<fptype> &x, const std::vector<fpcomplex> &y) -> std::vector<fpcomplex> {
    if(x.size() != y.size()) // Must be a valid pointer
        throw GeneralError("x and y must have the same diminsions!");

    int i, k;
    unsigned int n = x.size();
    std::vector<fpcomplex> u(n);
    std::vector<fpcomplex> y2(n);

    fptype sig, p, qn, un;
    fpcomplex yp1 = 2. * (y[1] - y[0]) / (x[1] - x[0]);
    fpcomplex ypn = 2. * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

    /* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/
    if(yp1.real() > 0.99e30) {
        y2[0].real(0.);
        u[0].real(0.);
    } else {
        y2[0].real(-0.5);
        u[0].real(3.0 / (x[1] - x[0]) * ((y[1].real() - y[0].real()) / (x[1] - x[0]) - yp1.real()));
    }
    if(yp1.imag() > 0.99e30) {
        y2[0].imag(0.);
        u[0].imag(0.);
    } else {
        y2[0].imag(-0.5);
        u[0].imag(3.0 / (x[1] - x[0]) * ((y[1].imag() - y[0].imag()) / (x[1] - x[0]) - yp1.imag()));
    }

    /* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the
     * decomposed factors*/

    for(i = 1; i < n - 1; i++) {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p   = sig * y2[i - 1].real() + 2.0;
        y2[i].real((sig - 1.0) / p);
        u[i].real((y[i + 1].real() - y[i].real()) / (x[i + 1] - x[i])
                  - (y[i].real() - y[i - 1].real()) / (x[i] - x[i - 1]));
        u[i].real((6.0 * u[i].real() / (x[i + 1] - x[i - 1]) - sig * u[i - 1].real()) / p);
        p = sig * y2[i - 1].imag() + 2.0;
        y2[i].imag((sig - 1.0) / p);
        u[i].imag((y[i + 1].imag() - y[i].imag()) / (x[i + 1] - x[i])
                  - (y[i].imag() - y[i - 1].imag()) / (x[i] - x[i - 1]));
        u[i].imag((6.0 * u[i].imag() / (x[i + 1] - x[i - 1]) - sig * u[i - 1].imag()) / p);
    }

    /* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

    if(ypn.real() > 0.99e30) {
        qn = 0.;
        un = 0.;
    } else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn.real() - (y[n - 1].real() - y[n - 2].real()) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1].real((un - qn * u[n - 2].real()) / (qn * y2[n - 2].real() + 1.0));
    if(ypn.imag() > 0.99e30) {
        qn = 0.;
        un = 0.;
    } else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn.imag() - (y[n - 1].imag() - y[n - 2].imag()) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1].imag((un - qn * u[n - 2].imag()) / (qn * y2[n - 2].imag() + 1.0));

    /* This is the backsubstitution loop of the tridiagonal algorithm */

    for(k = n - 2; k >= 0; k--) {
        y2[k].real(y2[k].real() * y2[k + 1].real() + u[k].real());
        y2[k].imag(y2[k].imag() * y2[k + 1].imag() + u[k].imag());
    }

    return y2;
}

__device__ fptype linear_interpolation(const fptype &s, const fptype *s0, const fptype *val, const size_t &nKnobs){

    //find the index
    int KhiAB = 0;

    while(KhiAB < nKnobs) {
        if(s < s0[KhiAB])
            break;
        KhiAB++;
    }

    if(KhiAB <= 0 || KhiAB == nKnobs)
        return 0.;

    int KloAB = KhiAB - 1;

    unsigned int twokloAB      = KloAB + KloAB;
    unsigned int twokhiAB      = KhiAB + KhiAB;
    fptype pwa_coefs_prime_real_kloAB = Deriatives[twokloAB];
    fptype pwa_coefs_prime_real_khiAB = Deriatives[twokhiAB];

    auto dmKK = s0[KhiAB] - s0[KloAB];
    auto aa = (s0[KhiAB] - s)/dmKK;
    auto bb = 1.0 - aa;
    auto aa3  = aa * aa * aa;
    auto bb3  = bb * bb * bb;

    //cubic interpolation
    auto coef_loAB = aa*val[KloAB] + (aa3 - aa)*pwa_coefs_prime_real_kloAB* (dmKK * dmKK) / 6.0;
    auto coef_hiAB = bb*val[KhiAB] + (bb3 - bb)*pwa_coefs_prime_real_khiAB* (dmKK * dmKK) / 6.0;
    coef_loAB+=coef_hiAB;
    return coef_loAB;
        // auto coef_loAB = val[KloAB];
    // auto coef_hiAB = val[KhiAB];
    // return aa*coef_loAB + bb*coef_hiAB; 
}

__device__ fpcomplex f0980_RBW(fptype s, fpcomplex scale){
        fptype resmass = 0.965;
        fptype gamma0 = 0.165;
        fptype resmass2 = resmass*resmass;

        fptype A = (resmass2 - s);
        fptype B = gamma0*resmass;
        fptype C = 1.0 / (A*A + B*B);

        fpcomplex ret(A * C, B * C);

        ret*= B;

        return ret*scale;
}

__device__ auto f0980_Flatte(fptype s, fpcomplex scale) -> fpcomplex {
    // indices[1] is unused constant index, for consistency with other function types.
    fptype resmass            = 0.984;
    fptype g1                 = 0.165;
    fptype g2                 = 4.21*g1;
    unsigned int particle     = 0;

    fpcomplex ret(0., 0.);

    if(resmass<0.){
        resmass *= -1.;
    }

    if(g1<0.){
        g1 *= -1.;
    }

    if(g2<0.){
        g2 *= -1.;
    }

    if(resmass<1.e-10 || g1<1.e-10 || g2<1.e-10){
        return ret;
    }

    fptype resmass2 = POW2(resmass);


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
    
    fptype rho1(0.0), rho2(0.0);
        
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
       
    
    return ret*scale ;
}

__device__ fpcomplex pelaez_pipi2pipi_interpolation(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) {
    fpcomplex ret(0.,0.);

    fptype akk_mag = pc.getConstant(0);
    fptype akk_phs = pc.getConstant(1);
    auto akk = thrust::polar(akk_mag,akk_phs);

    fptype apipi_mag = pc.getConstant(2);
    fptype apipi_phs = pc.getConstant(3);
    auto apipi = thrust::polar(apipi_mag,apipi_phs);

    fptype f0Scale_mag = pc.getConstant(4);
    fptype f0Scale_phs = pc.getConstant(5);
    auto f0Scale = thrust::polar(f0Scale_mag,f0Scale_phs);

    unsigned int charge_pos = pc.getConstant(6);
    unsigned int cyclic_index        = pc.getConstant(7);
    unsigned int symmDP       = pc.getConstant(8);
    unsigned int phi_nKnobs = pc.getConstant(9);

    //fpcomplex Lambda_d(-0.21874,2.51e-5); // Lambda_q = V_{cq}^*V_{uq} (this the D0bar CKM matrix element)
    auto Lambda_d = thrust::polar(0.21874000,3.1414779);
    // auto Lambda_s = thrust::polar(0.21890000,5.9387848e-06); 
    fpcomplex Lambda_s(0.21890,0.); 
    
    fptype sAB;
    fptype sAC;
    
    fpcomplex im(0.,1.);

    switch(cyclic_index) {
        case PAIR_12:
            sAB = m12;
            sAC = m13;
            break;
        case PAIR_13:
            sAB = m13;
            sAC = m23;
            break;
        case PAIR_23:
            sAB = m23;
            sAC = m13;
            break;
    }

    const fptype mk  = 0.496;
    // const fptype mpi  = 0.13957018; 
            

    #pragma unroll
    for(size_t i = 0; i < 1 + symmDP; i++) {

              
                if(charge_pos){
                    fptype eta_13 = Inela(sAB);
                    fptype arg_13 = argument(sAB)*M_PI/180.;
                    fpcomplex amp_pipi_13 =  eta_13*Lambda_d*thrust::polar(1.,2*arg_13)*apipi; 
                    
                    fptype eta_23 = Inela(sAC);
                    fptype arg_23 = argument(sAC)*M_PI/180.;
                    fpcomplex amp_pipi_23 =  eta_23*Lambda_d*thrust::polar(1.,2*arg_23)*apipi; 


                    fpcomplex amp_kk_13(0.,0.);    
                    fpcomplex amp_kk_23(0.,0.);  

                    if(sAB> 4.*mk*mk){

                        auto phase_KK_pipi = linear_interpolation(sAB, phi_s, phi_real, phi_nKnobs);
                        amp_kk_13 = im*sqrt(1.-eta_13*eta_13)*thrust::polar(1.,phase_KK_pipi)*Lambda_s*akk; 

                    }else if(sAC> 4.*mk*mk){

                        auto phase_KK_pipi = linear_interpolation(sAC, phi_s, phi_real, phi_nKnobs);
                        amp_kk_13 = im*sqrt(1.-eta_23*eta_23)*thrust::polar(1.,phase_KK_pipi)*Lambda_s*akk;

                    }else{
                        amp_kk_13 = fpcomplex(0.,0.);
                    }

                    //amp_pipi_13 *= f0980_Flatte(sAB,f0Scale);

                    fpcomplex amp = (amp_pipi_13+amp_kk_13); 

                    ret+=amp;

                   
                }else{

                    fptype eta_13 = Inela(sAB);
                    fptype arg_13 = argument(sAB)*M_PI/180.;
                    fpcomplex amp_pipi_13 =  eta_13*thrust::conj(Lambda_d)*apipi*thrust::polar(1.,2*arg_13); //amplitude pipi
                    
                    fptype eta_23 = Inela(sAC);
                    fptype arg_23 = argument(sAC)*M_PI/180.;
                    fpcomplex amp_pipi_23 =  eta_23*thrust::conj(Lambda_d)*apipi*thrust::polar(1.,2*arg_23); //amplitude pipi


                    fpcomplex amp_kk_13(0.,0.);    
                    fpcomplex amp_kk_23(0.,0.);  

                    if(sAB> 4.*mk*mk){

                        auto phase_KK_pipi = linear_interpolation(sAB, phi_s, phi_real, phi_nKnobs);
                        amp_kk_13 = im*sqrt(1.-eta_13*eta_13)*thrust::polar(1.,phase_KK_pipi)*thrust::conj(Lambda_s)*akk;

                    }else if(sAC> 4.*mk*mk){

                        auto phase_KK_pipi = linear_interpolation(sAC, phi_s, phi_real, phi_nKnobs);
                        amp_kk_13 = im*sqrt(1.-eta_23*eta_23)*thrust::polar(1.,phase_KK_pipi)*thrust::conj(Lambda_s)*akk;

                    }else{
                        amp_kk_13 = fpcomplex(0.,0.);
                    }

                    //amp_pipi_13 *= f0980_Flatte(sAB,f0Scale);

                    fpcomplex amp = (amp_pipi_13+amp_kk_13); 

                    ret+=amp;

                }


                fptype swpmass = sAB;
                sAB            = sAC;
                sAC            = swpmass;
    }
   
    
   
   return ret;
   
}

__device__ resonance_function_ptr ptr_to_AMP_PIPI = pelaez_pipi2pipi_interpolation;

namespace Resonances {
    ScatteringAmp::ScatteringAmp(std::string name,
                            Variable ar,
                            Variable ai,
                            std::complex<fptype> akk,
                            std::complex<fptype> apipi,
                            std::complex<fptype> f0scale,
                            std::vector<std::pair<fptype,std::complex<fptype>>> _phi00, // \delta_{KK}
                            unsigned int charge_pos,
                            unsigned int cyc,
                            bool symmDP)
            : ResonancePdf("ScatteringAmp",name, ar, ai){
        
            std::vector<unsigned int> pindices;

            //flags and constants

            registerConstant(abs(akk));
            registerConstant(arg(akk));

            registerConstant(abs(apipi));
            registerConstant(arg(apipi));

            registerConstant(abs(f0scale));
            registerConstant(arg(f0scale));

            registerConstant(charge_pos);
            registerConstant(cyc);
            registerConstant((unsigned int)symmDP);
            unsigned int phi_nKnobs = _phi00.size();
            printf("phi_nKnobs %d \n",phi_nKnobs);
            registerConstant(phi_nKnobs);

            registerFunction("ptr_to_AMP_PIPI", ptr_to_AMP_PIPI);

            //load deltaKK
            std::vector<fptype> host_phi_s;
            std::vector<fptype> host_phi_real;
            std::vector<fptype> host_phi_imag;
            std::vector<fpcomplex> host_phi;
            for(auto &v : _phi00){
                host_phi_s.push_back(v.first);
                host_phi_real.push_back(v.second.real());
                host_phi_imag.push_back(v.second.imag());
                host_phi.push_back(fpcomplex(v.second.real(),v.second.imag()));
            }
            MEMCPY_TO_SYMBOL(phi_s, host_phi_s.data(), host_phi_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);
            MEMCPY_TO_SYMBOL(phi_real, host_phi_real.data(), host_phi_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);
            MEMCPY_TO_SYMBOL(phi_imag, host_phi_imag.data(), host_phi_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);

            std::vector<fptype> y2_flat = make_flatten(calc_complex_derivative(host_phi_s, host_phi));
            MEMCPY_TO_SYMBOL(Deriatives, y2_flat.data(), 2 * phi_nKnobs * sizeof(fptype), 0, cudaMemcpyHostToDevice);

           
    }

}

}
