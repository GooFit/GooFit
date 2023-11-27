#include <goofit/PDFs/physics/resonances/ScatteringAmp.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {
__device__ fptype Deriatives[2 * 100];
__device__ fptype phi00_s[2*100];
__device__ fptype phi00_real[2*100];
__device__ fptype phi00_imag[2*100];

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

__device__ fpcomplex pelaez_pipi2pipi_interpolation(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) {
    fpcomplex ret(0.,0.);

    fptype akk = pc.getParameter(0);
    fptype apipi = pc.getParameter(1);
    unsigned int charge_pos = pc.getConstant(0);
    unsigned int cyclic_index        = pc.getConstant(1);
    unsigned int symmDP       = pc.getConstant(2);
    unsigned int phi00_nKnobs = pc.getConstant(3);

    fpcomplex Lambda_d(-0.21874,-2.51e-5); // Lambda_q = V_{cq}V_{uq}^* (this the D0bar CKM matrix element)
    fpcomplex Lambda_s(0.21890,-0.13e-5);
    
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

    fptype mAB = sqrt(sAB), mAC = sqrt(sAC);

    #pragma unroll
    for(size_t i = 0; i < 1 + symmDP; i++) {

                fptype phase_KK_pipi = 0.; //phase \delta_{KK}
                fpcomplex amp_kk = 0.; //amplitude KK
                
                if(sAB>phi00_s[0])
                    phase_KK_pipi = linear_interpolation(sAB, phi00_s, phi00_real, phi00_nKnobs);

                if(charge_pos){
                    fpcomplex eta(Inela(sAB),0.);
                    fpcomplex amp_pipi = eta*thrust::conj(Lambda_d)*thrust::exp(2.*im*argument(sAB)*M_PI/180.)*apipi; //amplitude pipi
                    if(sAB>phi00_s[0])
                        amp_kk = im*thrust::sqrt(1.-eta*eta)*thrust::exp(im*(phase_KK_pipi))*thrust::conj(Lambda_s)*akk; 
                    fpcomplex amp = amp_pipi + amp_kk; 
                    ret+=amp;
                }else{
                    fpcomplex eta(Inela(sAB),0.);
                    fpcomplex amp_pipi = eta*Lambda_d*thrust::exp(2.*im*argument(sAB)*M_PI/180.)*apipi; //amplitude pipi
                    if(sAB>phi00_s[0])
                        amp_kk = im*thrust::sqrt(1.-eta*eta)*thrust::exp(im*(phase_KK_pipi))*Lambda_s*akk; 
                    fpcomplex amp = amp_pipi + amp_kk; 
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
                            Variable akk,
                            Variable apipi,
                            std::vector<std::pair<fptype,fpcomplex>> &_phi00, // \delta_{KK}
                            bool charge_pos,
                            unsigned int cyc,
                            bool symmDP)
            : ResonancePdf("ScatteringAmp",name, ar, ai) {
        
            std::vector<unsigned int> pindices;

            //flags and constants
            registerParameter(akk);
            registerParameter(apipi);
            registerConstant(charge_pos);
            registerConstant(cyc);
            registerConstant((unsigned int)symmDP);
            unsigned int phi00_nKnobs = _phi00.size();
            printf("phi00_nKnobs %d \n",phi00_nKnobs);
            registerConstant(phi00_nKnobs);

            //load deltaKK
            std::vector<fptype> host_phi00_s;
            std::vector<fptype> host_phi00_real;
            std::vector<fptype> host_phi00_imag;
            std::vector<fpcomplex> host_phi;
            for(auto &v : _phi00){
                host_phi00_s.push_back(v.first);
                host_phi00_real.push_back(v.second.real());
                host_phi00_imag.push_back(v.second.imag());
                host_phi.push_back(fpcomplex(v.second.real(),v.second.imag()));
            }
            MEMCPY_TO_SYMBOL(phi00_s, host_phi00_s.data(), host_phi00_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);
            MEMCPY_TO_SYMBOL(phi00_real, host_phi00_real.data(), host_phi00_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);
            MEMCPY_TO_SYMBOL(phi00_imag, host_phi00_imag.data(), host_phi00_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);

            std::vector<fptype> y2_flat = make_flatten(calc_complex_derivative(host_phi00_s, host_phi));
            MEMCPY_TO_SYMBOL(Deriatives, y2_flat.data(), 2 * phi00_nKnobs * sizeof(fptype), 0, cudaMemcpyHostToDevice);

            registerFunction("ptr_to_AMP_PIPI", ptr_to_AMP_PIPI);
    }
}
}
