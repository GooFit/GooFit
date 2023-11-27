#include <goofit/PDFs/physics/resonances/ScatteringAmp.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ fptype phi00_s[2*100];
__device__ fptype phi00_real[2*100];
__device__ fptype phi00_imag[2*100];

///////////Amp PIPI->PIPI equation 7b (PHYSICAL REVIEW LETTERS 131, 051802 (2023))

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

    auto dmKK = s0[KhiAB] - s0[KloAB];
    auto aa = (s0[KhiAB] - s)/dmKK;
    auto bb = 1.0 - aa;

    auto coef_loAB = val[KloAB];
    auto coef_hiAB = val[KhiAB];

    return aa*coef_loAB + bb*coef_hiAB; 

}

__device__ fpcomplex pelaez_pipi2pipi_interpolation(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    fpcomplex ret(0.,0.);

    fptype akk = pc.getParameter(0);
    fptype apipi = pc.getParameter(1);
    unsigned int charge_pos = pc.getConstant(0);
    unsigned int cyclic_index        = pc.getConstant(1);
    unsigned int symmDP       = pc.getConstant(2);
    unsigned int phi00_nKnobs = pc.getConstant(3);

    fpcomplex Lambda_d(-0.21874,-2.51e-5);
    fpcomplex Lambda_s(0.21890,-0.13e-5);

    fpcomplex Vcd(-0.2,0.);
    fpcomplex Vud(1.,0.);
    fpcomplex Vcs(1.,0.);
    fpcomplex Vus(1.,0.);
    
    fptype sAB = m12, sAC = m13;
    fptype mAB = sqrt(m12), mAC = sqrt(m13);

    fpcomplex im(0.,1.);

    switch(cyclic_index) {
        case PAIR_13:
            sAB = m13;
            sAC = m12;
            break;
        case PAIR_23:
            sAB = m23;
            sAC = m12;
            break;
    }

    //printf("Phi(0.8)=%f \n",Phi(0.8));
    // fpcomplex tl = t00(0.8);
    // printf("t00=%f + %f i \n",tl.real(),tl.imag());

    //#pragma unroll
    for(size_t i = 0; i < 1 + symmDP; i++) {

                fptype phase_KK_pipi = 0.; //phase \delta_{KK}
                fpcomplex amp_kk = 0.; //amplitude KK
                
                if(sAB>1.)
                    phase_KK_pipi = linear_interpolation(sAB, phi00_s, phi00_real, phi00_nKnobs);

                if(charge_pos){
                    fpcomplex eta(Inela(sAB),0.);
                    fpcomplex amp_pipi = eta*thrust::conj(Lambda_d)*thrust::exp(2.*im*argument(sAB)*M_PI/180.)*apipi; //amplitude pipi
                    if(sAB>1.02)
                        amp_kk = im*thrust::sqrt(1.-eta*eta)*thrust::exp(im*(phase_KK_pipi))*thrust::conj(Lambda_s)*akk; 
                    fpcomplex amp = amp_pipi + amp_kk; 
                    ret+=amp;
                }else{
                    fpcomplex eta(Inela(sAB),0.);
                    fpcomplex amp_pipi = eta*Lambda_d*thrust::exp(-2.*im*argument(sAB)*M_PI/180.)*apipi; //amplitude pipi
                    if(sAB>1.02)
                        amp_kk = -im*thrust::sqrt(1.-eta*eta)*thrust::exp(-im*(phase_KK_pipi))*Lambda_s*akk; 
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
            for(auto &v : _phi00){
                host_phi00_s.push_back(v.first);
                host_phi00_real.push_back(v.second.real());
                host_phi00_imag.push_back(v.second.imag());
            }
            MEMCPY_TO_SYMBOL(phi00_s, host_phi00_s.data(), host_phi00_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);
            MEMCPY_TO_SYMBOL(phi00_real, host_phi00_real.data(), host_phi00_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);
            MEMCPY_TO_SYMBOL(phi00_imag, host_phi00_imag.data(), host_phi00_s.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);

            registerFunction("ptr_to_AMP_PIPI", ptr_to_AMP_PIPI);
    }
}
}
