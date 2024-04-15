#include <goofit/PDFs/physics/resonances/Rescattering2.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ auto pn( const fptype x_, const fptype n) -> fptype  {

    if (n==0) return 1.0;
	if (n==1) return x_;
	return 2*x_*pn(x_,n-1) - pn(x_,n-2);

}

__device__ auto x(const fptype sqr_t, const fptype *sqr_tmin, const fptype *sqr_tmax, const int i) -> fptype  
{       
	return 2.0*(sqr_t-sqr_tmin[i])/(sqr_tmax[i]-sqr_tmin[i]) -1.0;
}

__device__ auto phi00(const fptype sqr_t, const fptype *sqr_tmin, const fptype *sqr_tmax,  const int i, const fptype* coefs) -> fptype  
{       

    auto B0_ = coefs[0];
    auto B1_ = coefs[1];
    auto B2_ = coefs[2];
    auto B3_ = coefs[3];
    auto C0_ = coefs[4];
    auto C1_ = coefs[5];
    auto C2_ = coefs[6];
    auto C3_ = coefs[7];
    auto C4_ = coefs[8];
    auto C5_ = coefs[9];


	fptype x_t = x(sqr_t, sqr_tmin, sqr_tmax, i);

	if (i==1) return                 B0_*pn(x_t,0)+ 
				B1_*pn(x_t,1)+
				B2_*pn(x_t,2)+
				B3_*pn(x_t,3);

	if (i==2) return                 C0_*pn(x_t,0)+ 
				C1_*pn(x_t,1)+
				C2_*pn(x_t,2)+
				C3_*pn(x_t,3)+
				C4_*pn(x_t,4)+
				C5_*pn(x_t,5);

	return 0;
}

__device__ auto g00(const fptype sqr_t, const fptype *sqr_tmin,const fptype *sqr_tmax, const int i, const fptype* coefs) -> fptype  
{       
	fptype x_t = x(sqr_t, sqr_tmin,sqr_tmax, i);

    auto D0_ = coefs[0];
    auto D1_ = coefs[1];
    auto D2_ = coefs[2];
    auto D3_ = coefs[3];
    auto F0_ = coefs[4];
    auto F1_ = coefs[5];
    auto F2_ = coefs[6];
    auto F3_ = coefs[7];
    auto F4_ = coefs[8];


	if (i==1) return        D0_*pn(x_t,0)+ 
				D1_*pn(x_t,1)+
				D2_*pn(x_t,2)+
				D3_*pn(x_t,3);


	if (i==2) return                 F0_*pn(x_t,0)+ 
				F1_*pn(x_t,1)+
				F2_*pn(x_t,2)+
				F3_*pn(x_t,3)+
				F4_*pn(x_t,4);

	return 0;
}

__device__ auto qh(const fptype s, const fptype m) -> fptype 
{ //p9 
  if ( s < 4. * m * m){
  return 0.0;
  }

  if (s >= 4. * m * m){
    fptype q = (sqrt(s-4*m*m) /2. );
  return q;
  }

  return 0;

}
 

template <int I>
__device__ auto func_Rescattering2(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex {
  
    unsigned int cyclic_index = pc.getConstant(0);

    fptype B1_ = pc.getParameter(0);
    fptype B2_ = pc.getParameter(1);
    fptype B3_ = pc.getParameter(2);
    fptype C1_ = pc.getParameter(3);
    fptype C2_ = pc.getParameter(4);
    fptype C3_ = pc.getParameter(5);
    fptype C4_ = pc.getParameter(6);
    fptype C5_ = pc.getParameter(7);
    fptype D0_ = pc.getParameter(8);
    fptype D1_ = pc.getParameter(9);
    fptype D2_ = pc.getParameter(10);
    fptype D3_ = pc.getParameter(11);
    fptype F1_ = pc.getParameter(12);
    fptype F2_ = pc.getParameter(13);
    fptype F3_ = pc.getParameter(14);
    fptype F4_ = pc.getParameter(15);

    fptype B0_ = 0.;
    fptype C0_ = 295.39;
    fptype F0_ = 0.151;
 
    fptype s  = 0.0;
    fptype mass = 0.0;

    fpcomplex ret(0.,0.);

    const fptype k_MASS = 0.493677;
    const fptype sqr_tmin[3] = {0.,2.*k_MASS,sqrt(2.)};
    const fptype sqr_tmax[3] = {0.,sqrt(2.),2.};

#pragma unroll
    for(size_t j = 0; j < I; j++) {
        
        if(PAIR_12 == cyclic_index){
            s    = m12 ;
        }

        if(PAIR_13 == cyclic_index){
            s    = m13 ;
        }

         if(PAIR_23 == cyclic_index){
            s    = m23 ;
        }
       
        mass = sqrt(s);
  

        B0_ = 226.5 + B1_ - B2_ + B3_;

        const fptype coefs_phi00[10] = {B0_,B1_,B2_,B3_,C0_,C1_,C2_,C3_,C4_,C5_};

        //C0_ = phi00(sqr_tmax[1],sqr_tmin, sqr_tmax, 1, coefs_phi00) + C1_ - C2_ + C3_ - C4_ + C5_;

        const fptype coefs_g00[9] = {D0_,D1_,D2_,D3_,F0_,F1_,F2_,F3_,F4_};

        //F0_ = g00(sqr_tmax[1], sqr_tmin, sqr_tmax, 1, coefs_g00) + F1_ - F2_ + F3_ - F4_;

        // printf("g00(sqr_tmax[1],1)= %f \n",g00(sqr_tmax[1], sqr_tmin, sqr_tmax, 1, coefs_g00));
        // printf("g00(sqr_tmin[1],1)= %f \n",g00(sqr_tmin[1], sqr_tmin, sqr_tmax, 1, coefs_g00));
        // printf("g00(sqr_tmax[2],2)= %f \n",g00(sqr_tmax[2], sqr_tmin, sqr_tmax, 2, coefs_g00));
        // printf("ph00(sqr_tmax[1],1)= %f \n",phi00(sqr_tmax[1], sqr_tmin, sqr_tmax, 1, coefs_phi00));
        // printf("ph00(sqr_tmin[1],1)= %f \n",phi00(sqr_tmin[1], sqr_tmin, sqr_tmax, 1, coefs_phi00));
        // printf("ph00(sqr_tmax[2],2)= %f \n",phi00(sqr_tmax[2], sqr_tmin, sqr_tmax, 2, coefs_phi00));


        int i=0 ; 
        if (mass < sqr_tmax[1]) i = 1;
        if (mass > sqr_tmax[1] && mass < sqr_tmax[2]) i = 2;
        if (i == 0) {
            return fpcomplex(0,0);
        }
        fptype NR1_s = 1.0/(1.0+s);
        fptype mag   = g00(mass, sqr_tmin, sqr_tmax, i, coefs_g00);
        if (mass > sqr_tmin[2]) mag=0.0;
        fptype phase = phi00(mass, sqr_tmin, sqr_tmax, i, coefs_phi00)*M_PI/180.0;

        ret+= fpcomplex(mag*cos(phase)*NR1_s, mag*sin(phase)*NR1_s);

         if(I != 0) {
            fptype swpmass = m13;
            m13            = m23;
            m23            = swpmass;
        }
    }

    
    pc.incrementIndex(1, 16, 1, 0, 1);
  
    return ret;
}

__device__ resonance_function_ptr ptr_to_Rescattering2     = func_Rescattering2<1>;
__device__ resonance_function_ptr ptr_to_Rescattering2_Sym = func_Rescattering2<2>;

namespace Resonances {

Rescattering2::Rescattering2(std::string name,
         Variable ar,
         Variable ai,
         std::vector<Variable> coefs,
         unsigned int cyc,
         bool sym)
    : ResonancePdf("Rescattering2", name, ar, ai) {

    if(coefs.size()!=16)
        GOOFIT_ERROR("ERROR: Rescattering2 expects 16 coefs");

    for(int i=0; i<coefs.size();i++)
        registerParameter(coefs.at(i));

    registerConstant(cyc);



    if(sym)
        registerFunction("ptr_to_Rescattering2_Sym", ptr_to_Rescattering2_Sym);
    else
        registerFunction("ptr_to_Rescattering2", ptr_to_Rescattering2);
}

} // namespace Resonances

} // namespace GooFit
