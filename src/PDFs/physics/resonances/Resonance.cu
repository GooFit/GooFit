#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

    __device__ auto h(const fptype &m,const fptype &q)->fptype{
        auto const mpi = 0.13957018;
        return (2.*q/(M_PI*m))*log( (m+2.*q)/(2.*mpi) );
    }
    
    __device__ auto h_prime(const fptype &m0,const fptype &q0)->fptype{
        return (h(m0,q0)*( (1./(8.*q0*q0)) - (1./(2.*m0*m0) ) + (1./(2.*M_PI*m0*m0))));
    }
    
    __device__ auto d(const fptype &m0,const fptype &q0)->fptype{
        auto const mpi = 0.13957018;
        return ((3.*POW2(mpi)/(M_PI*POW2(q0)))*log( (m0+2.*q0)/(2.*mpi)) + (m0/(2.*M_PI*q0)) - (POW2(mpi)*m0/(M_PI*POW2(q0)*q0)));
    }
    
    __device__ auto f(const fptype &m, const fptype &m0,const fptype &width , const fptype &q, const fptype &q0)->fptype{
        return width*(POW2(m0)/(POW2(q0)*q0))*(POW2(q)*(h(m,q)-h(m0,q0)) + (POW2(m0)-POW2(m))*q0*q0*h_prime(m0,q0));
    }


__device__ auto calc_q(fptype s12, fptype m1, fptype m2)-> fptype{

    fptype EiCmsij = (s12 - m1*m1 + m2*m2)/(2.*sqrt(s12));

    auto arg = (EiCmsij*EiCmsij - m2*m2);

    if(arg<0.0){
        return 0.0;
    }

    double q = sqrt(arg);

    return q;
}


__device__ auto DaugDecayMomResFrame(fptype rMassSq, fptype d1m, fptype d2m) -> fptype {
  // Decay momentum of either daughter in the resonance rest frame
	// when resonance mass = rest-mass value, m_0 (PDG value)

    fptype term1 = rMassSq - (d1m+d2m)*(d1m+d2m);
    fptype term2 = rMassSq - (d1m-d2m)*(d1m-d2m);
    fptype term12 = term1*term2;
    fptype q      = 1.0;

    //printf("%f \t %f  \t %f  \t %f  \t %f \n ",term1,term2,rMassSq,d1m,d2m); 
    if(term12> 0.0){
        q = sqrt(term12)/(2.0*sqrt(rMassSq));
    }else{
        GOOFIT_TRACE("term12 < zero!");
        q = 0.0;
    }


    return q;
}



__device__ auto BachMomResFrame(fptype M, fptype rMassSq, fptype mBach) -> fptype {
    // Momentum of the bachelor particle in the resonance rest frame
	// when resonance mass = rest-mass value, m_0 (PDG value)
  
      fptype eBach = (M*M - rMassSq - mBach*mBach)/(2.0*sqrt(rMassSq));
      fptype termBach = eBach*eBach - mBach*mBach;
      fptype p      = 1.0;

        if ( eBach<0.0 || termBach<0.0 ) {
            p = 0.0;
            GOOFIT_TRACE("eBach<0.0 || termBach<0.0");
        } else {
            p = sqrt( termBach );
        }

        return p;
  }

  __device__ auto BachMomParentFrame(fptype M, fptype mBach, fptype rMassSq) -> fptype {
   	// Momentum of the bachelor particle in the parent rest frame
	// when resonance mass = rest-mass value, m_0 (PDG value)
  
      fptype eStarBach = (M*M + mBach*mBach - rMassSq)/(2.0*M);
      fptype termStarBach = eStarBach*eStarBach - mBach*mBach;
      fptype pstar      = 1.0;

        if ( eStarBach<0.0 || termStarBach<0.0 ) {
            pstar = 0.0;
            GOOFIT_TRACE("eStarBach<0.0 || termStarBach<0.0");
        } else {
            pstar = sqrt( termStarBach );
        }

        return pstar;
  }



  __device__ auto BlattWeisskopfPrime(fptype z, unsigned int spin)-> fptype {

        fptype ret = 1.;

        if(spin==0)
            return 1.0;
            
        switch (spin){
            case 1:
                ret = 1.0/sqrt(1.0 + z*z);
                break;
            case 2:
                ret = 1./sqrt(z*z*z*z + 3.0*z*z + 9.0);
                break;

            case 3:
                ret = 1./sqrt(z*z*z*z*z*z + 6.0*z*z*z*z + 45.0*z*z + 255.0);
                break;
        }

        return ret;

  }

  
__device__ auto cFromM(
    fptype motherMass,
    fptype daug1Mass,
    fptype daug2Mass,
    fptype daug3Mass,
    fptype m13,
    fptype m23,
    fptype m12,
    unsigned int cyclic_index) -> fptype {

        fptype  _mA  = 0.;
        fptype  _mB  = 0.;
        fptype  _mC  = 0.;
        fptype  _mAC = 0.;
        fptype  _mAB = 0.;

        if(PAIR_12 == cyclic_index){
            _mAB = m12;
            _mAC = m13;
            _mA = daug1Mass;
            _mB = daug2Mass;
            _mC = daug3Mass;
        }

        if(PAIR_13 == cyclic_index){
            _mAB = m13;
            _mAC = m23;
            _mA = daug3Mass;
            _mB = daug1Mass;
            _mC = daug2Mass;
        }

         if(PAIR_23 == cyclic_index){
            _mAB = m23;
            _mAC = m12;
            _mA = daug2Mass;
            _mB = daug3Mass;
            _mC = daug1Mass;
        }
    

        fptype EACmsAB = (_mAB - _mB*_mB + _mA*_mA)/(2.0*sqrt(_mAB));
        fptype ECCmsAB = (motherMass*motherMass - _mAB - _mC*_mC)/(2.0*sqrt(_mAB));

        if(EACmsAB<_mA){
            GOOFIT_TRACE("EACmsAB<_mA");
            return 0.0;
        }


        if(ECCmsAB<_mC){
            GOOFIT_TRACE("ECCmsAB<_mC");
            return 0.0;
        }


        fptype qA_ = sqrt(EACmsAB*EACmsAB - _mA*_mA);
        fptype qC_ = sqrt(ECCmsAB*ECCmsAB - _mC*_mC);

        fptype cosHel = -(_mAC - _mA*_mA - _mC*_mC - 2.0*EACmsAB*ECCmsAB)/(2.0*qA_*qC_);

        if(cosHel > 1.0){
            cosHel = 1.0;
        }else if(cosHel<-1.0){
            cosHel = -1.0;
        }

        if(cyclic_index==PAIR_12 || cyclic_index==PAIR_13 )
            cosHel *= -1.0;

        return cosHel;

}

__device__ auto calcLegendrePoly(fptype cosHel, unsigned int spin) -> fptype{
    fptype legPol = 1.0;

    if(spin==0)
        return 1.0;

    switch(spin){
        case 1:
            legPol = -2.0*cosHel;
            break;
        case 2:
            legPol = 4.0*(3.0*cosHel*cosHel - 1.0)/3.0;
            break;
        case 3:
            legPol = -8.0*(5.0*cosHel*cosHel*cosHel - 3.0*cosHel)/5.0;
            break;
    }

    return legPol;
}

__device__ auto calcZemachSpinFactor(fptype pProd, fptype legPol, unsigned int spin) -> fptype{
 

    if(spin==0)
        return 1.0;

    fptype spinFactor(pProd);

    for ( int i(1); i < spin; ++i ) {
        spinFactor *= pProd;
    }
    
    spinFactor *= legPol;
    
    return spinFactor;
}


  
  __device__ auto twoBodyCMmom(double rMassSq, fptype d1m, fptype d2m) -> fptype {
    // For A -> B + C, calculate momentum of B and C in rest frame of A.
    // PDG 38.16.

    fptype kin1 = 1 - POW2(d1m + d2m) / rMassSq;

    kin1 = kin1 >= 0 ? sqrt(kin1) : 1;

    fptype kin2 = 1 - POW2(d1m - d2m) / rMassSq;
    kin2        = kin2 >= 0 ? sqrt(kin2) : 1;

    return 0.5 * sqrt(rMassSq) * kin1 * kin2;
}

__device__ auto twoBodyCMMothermom(fptype rMassSq, fptype dm, fptype d3m) -> fptype {
    fptype kin1 = 1 - POW2(dm + d3m) / rMassSq;
    if(kin1 >= 0)
        kin1 = sqrt(kin1);
    else
        kin1 = 1;
    fptype kin2 = 1 - POW2(dm - d3m) / rMassSq;
    if(kin2 >= 0)
        kin2 = sqrt(kin2);
    else
        kin2 = 1;

    return 0.5 * rMassSq * kin1 * kin2 / dm;
}

__device__ auto dampingFactorSquare(const fptype &cmmom, const int &spin, const fptype &mRadius) -> fptype {
    fptype square = mRadius * mRadius * cmmom * cmmom;
    fptype dfsq   = 2 * square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    fptype dfsqres = 13 * square / pow(square - 3, 2) + 9 * square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}

__device__ auto dampingFactorSquareNorm(const fptype &cmmom, const int &spin, const fptype &mRadius) -> fptype {
    fptype square = mRadius * mRadius * cmmom * cmmom;
    fptype dfsq   = 1 + square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    fptype dfsqres = dfsq + 8 + 2 * square + square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}

__device__ auto spinFactor(unsigned int spin,
                           fptype motherMass,
                           fptype daug1Mass,
                           fptype daug2Mass,
                           fptype daug3Mass,
                           fptype m12,
                           fptype m13,
                           fptype m23,
                           unsigned int cyclic_index) -> fptype {
    auto ret = 1.;

    auto const _mA  = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
    auto const _mB  = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
    auto const _mC  = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
    auto const _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12));
    auto const _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m12 : m13));
    auto const _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

    if(1 == spin) {
        auto const massFactor = 1.0 / _mAB;
        ret = ((_mBC - _mAC) + (massFactor * (motherMass * motherMass - _mC * _mC) * (_mA * _mA - _mB * _mB)));
    }

    if(2 == spin) {
        auto const massFactor = 1.0 / _mAB;
        auto const a1
            = ((_mBC - _mAC) + (massFactor * (motherMass * motherMass - _mC * _mC) * (_mA * _mA - _mB * _mB)));
        auto const a2 = ((_mAB - (2 * motherMass * motherMass) - (2 * _mC * _mC))
                         + massFactor * POW2(motherMass * motherMass - _mC * _mC));
        auto const a3 = ((_mAB - (2 * _mA * _mA) - (2 * _mB * _mB)) + massFactor * POW2(_mA * _mA - _mB * _mB));

        ret = POW2(a1) - a2 * a3 / 3;
    }

    return ret;
}


//Defining constants 
__constant__ fptype EPS = 1.0e-10;
__constant__ fptype mk  = 0.496;
__constant__ fptype mp  = 0.13957018;

//defining arrays n==1

__constant__   int ni = 1;

__constant__  fptype b00 = 12.23965 ;
__constant__  fptype b01 = -0.9395; 
__constant__  fptype b02 = 15.8534; 
__constant__  fptype b03 = -5.7445; 
__constant__  fptype b04 = -22.4524; 
__constant__  fptype b05 = 6.92; 
__constant__  fptype b06 = 0.; 
__constant__  fptype b07 = 0.; 
__constant__  fptype k00 = 1.28845; 
__constant__  fptype k01 = -1.08423; 
__constant__  fptype k02 = 0.043397;
__constant__  fptype k03 = -0.067912;
__constant__  fptype k04 = 0.;
__constant__  fptype mf00 = 0.9963;
__constant__  fptype wf00= -0.0251;
__constant__  fptype d00= -5.3625;
__constant__  fptype d01 =0.;
__constant__  fptype d02 =0.;
__constant__  fptype e00= 10.34268;
__constant__  fptype e01= 0.;
__constant__  fptype e02= 0.;
__constant__  fptype z00= 0.136615;
__constant__  fptype sm = 1.96; 


__device__ fpcomplex _Vc (fptype c){ //make a variable fptype and constant became complex;
    return fpcomplex(c,0.);
}

__device__ fpcomplex sigma (fpcomplex s, fptype m){
    
    if(s.imag()>=0.){
        return thrust::sqrt(1. - (4.*m*m)/s);
    }else{
       return - thrust::sqrt(1. - (4.*m*m)/s);
    }
}

__device__ fpcomplex V(fpcomplex x){
     fpcomplex s0(4.*mk*mk,0.);
     return
     (thrust::sqrt(x)-thrust::sqrt(s0-x))/(thrust::sqrt(x)+thrust::sqrt(s0-x));
}

__device__ fpcomplex wp(fpcomplex x){
    const fpcomplex epss(0, 1.*pow(10,-18));
    return V(x + epss); 
}

__device__ fptype q(fptype s){
    return sqrt(s/4. - mp*mp);
}

__device__ fpcomplex tlow(fpcomplex s){
    fpcomplex c1((-(z00*z00)/2.),0.);
    fpcomplex c2(z00*z00/(mp),0.);
    fpcomplex phi00sig = thrust::sqrt(s)*(_Vc(mp*mp)/(c1+s)*(_Vc(b00) + c2/thrust::sqrt(s) + _Vc(b01)*wp(s) + _Vc(b02)*wp(s)*wp(s)  + _Vc(b03)*wp(s)*wp(s)*wp(s) + _Vc(b04)*wp(s)*wp(s)*wp(s)*wp(s) + _Vc(b05)*wp(s)*wp(s)*wp(s)*wp(s)*wp(s)))/_Vc(2.*q(s.real()));
   
    fpcomplex phi00= thrust::sqrt(s)*_Vc(mp*mp)*((c2/thrust::sqrt(s)) + _Vc(b00) + _Vc(b01)*wp(s) + _Vc(b02)*wp(s)*wp(s)  + _Vc(b03)*wp(s)*wp(s)*wp(s) + _Vc(b04)*wp(s)*wp(s)*wp(s)*wp(s) + _Vc(b05)*wp(s)*wp(s)*wp(s)*wp(s)*wp(s))/(_Vc(2.*q(s.real()))*(s+c1));    
    fpcomplex c3( 0., 1.);
    fpcomplex low =fpcomplex(1.,0.)/(sigma(s,mp)*(phi00 - c3));
    return low;
}

__device__ fpcomplex pn( fpcomplex x, fptype n) {
     if (n==0.) return fpcomplex(1.,0.);
     if (n==1.) return x;
    fpcomplex C(2.,0.);
    return C*x*pn(x,n-1.) -pn(x,n-2.);
}

__device__ fpcomplex fu(fpcomplex x){
     fpcomplex w1=_Vc(2.)*(thrust::sqrt(x) - _Vc(2.*mk))/_Vc(1.5-2.*mk) -fpcomplex(1.,0.) ;
     return (_Vc(k00) + _Vc(k01)*pn(w1,1.) + _Vc(k02)*pn(w1,2.) + _Vc(k03)*pn(w1,3.));
}
__device__ fpcomplex Jp(fpcomplex s,fptype m){
    fpcomplex c1(2.0/M_PI,0.);
    fpcomplex c2(1.0/M_PI,0.);
    fpcomplex a=sigma(s,m);
    fpcomplex alpha = (a- fpcomplex(1.,0.))/(a + fpcomplex(1.,0.));
    fpcomplex J = c1+ c2*a*thrust::log(alpha);
    return J;
}

__device__ fpcomplex tf0(fpcomplex s){

    const fptype sR = 0.991984;
    const fptype sI = - 0.0500143;
    const fpcomplex spol(sR,sI);

    fptype fR = fu(spol).real();
    fptype fI = fu(spol).imag();

    fptype Jpr= Jp(spol,mp).real() ; //-0.554509;
    fptype Jpi=  Jp(spol,mp).imag() ; // -0.941444;
    
    fptype Jkr= Jp(spol,mk).real() ;// 0.48917;//
    fptype Jki= Jp(spol,mk).imag() ;//-0.145356;//
    
    fptype sigR= sigma(spol,mp).real();
    fptype sigI= sigma(spol,mp).imag();

    fptype d = Jpi*sR + Jpr*sI + 2.*(sI*sigI - sR*sigR);
    fpcomplex G (- (fI*Jkr + fR*Jki + sI)/d,0.);
    fpcomplex M (( (fI*Jkr + fR*Jki)*(sI*(Jpi - 2.*sigR) - sR*(Jpr + 2.*sigI)) + (Jpi - 2.*sigR)*(sI*sI + sR*sR))/d - (fI*Jki - fR*Jkr),0.);
    
    fpcomplex jpp = Jp(s,mp);
    fpcomplex jpk = Jp(s,mk);
    fpcomplex Fu = fu(s);
    return (s*G)/(M - s - s*G*jpp- Fu*jpk);
}

//up to here all without interference (But Is different from Mathematica)
__device__ fpcomplex t00(fptype s){
    fpcomplex C1(0.,2.); // in the Mathematica code is defined as -sigma to all s real. ==> Key 1
    fpcomplex tl= tlow(s);
    fpcomplex tf = tf0(s);
    fpcomplex sig = sigma(s,mp);
  return tl + tf + C1*sig*tl*tf;
}

__device__ fpcomplex S00(fptype s){
    fpcomplex C(0.,2.); // in the Mathematica code is defined as -sigma to all s real. ==> Key 2
    fpcomplex sig = sigma(s,mp);
    fpcomplex t = t00(s);
    return fpcomplex(1.,0.) + C*sig*t;
}

//Inenastic is need to Energy>1.4=sm (eq 16 - Global P)
__device__ fptype wp2(fptype s){
    fptype qm = q(sm);
    fptype tf12 = sqrt(sm);
    fptype tf22 = 2.;
    return 2.*sqrt(s)/(tf22 - tf12) - (2.*tf12)/(tf22 - tf12) - 1.;
}

__device__ fptype Phim(fptype s){
    fptype fase= thrust::arg(S00(s)); //corrected
    return fase*(90./M_PI) + 360.; // checar esse angulo!
}

__device__ fptype derivaPhi(fptype s){ //deri, derip, deripp
    fpcomplex ret;
    // if(ni==1){
        ret = (96.3 - d00*((pn(wp2(s),2) - pn(wp2(s)-EPS,2))/EPS) )/(( pn(wp2(s),1) - pn(wp2(s)-EPS,1))/EPS );
        return ret.real();
    // }
    // if(ni==2){
    //     return (95.65 - d00*(( pn(wp2(s),2) - pn(wp2(s)-EPS,2))/EPS) )/((pn(wp2(s),1) - pn(wp2(s)-EPS,1))/EPS).real(); //s -> sm;
    // }
    // if(ni==3){
    //     return (94.2749 - d00*( (pn(wp2(s),2) - pn(wp2(s)-EPS,2))/EPS) - d01*(pn(wp2(s),3) - pn(wp2(s)-EPS,3)).real()/EPS - d02*(pn(wp2(s),4) - pn(wp2(s)-EPS,4)).real()/EPS)/((pn(wp2(s),1) - pn(wp2(s)-EPS,1)).real()/EPS);
    // }

    
}

__device__ fptype Phi(fptype s){
    // if(ni==1){
        return (Phim(sm) + derivaPhi(sm)*(pn(wp2(s),1) + 1.) + d00*(pn(wp2(s),2) - 1.)).real();
    // }
    // if(ni==2){
    //     return (Phim(sm) + derivaPhi(sm)*(pn(wp2(s),1) + 1.) +
    //     d00*(pn(wp2(s),2) - 1.)).real();
    // }
    // if(ni==3){
    //     return (Phim(sm) + derivaPhi(sm)*(pn(wp2(s),1) + 1.) +
    //     d00*(pn(wp2(s),2) - 1.) + d01*(pn(wp2(s),3) + 1.) + d02*(pn(wp2(s),4) - 1.)).real();
    // }
}
__device__ fptype argument(fptype s){
    fptype ImS00 =S00(s).imag();
    fptype ReS00=S00(s).real();
    fptype fase =(90.0/M_PI)*atan2(ImS00,ReS00);

    if(s<sm){
        if(s<0.81 && fase>0.){
            return fase ;
        } else if(s<1.21 || fase>0. ){
            return fase + 180.;
        } else{
            return fase + 360.;
        }
    }else{ 
        return Phi(s);
    }
}

__device__ fptype Inela(fptype s){
    fptype deri2 = 0.61996;
    fptype deri2p = 1.3401840137492784;
    fptype deri2pp = 0.8728282172774287;
    fptype S00abs=thrust::abs(S00(sm)); //sm=1.96 and s00abs=3.071756
    //printf("S00abs = %f \n",S00abs);
    fptype Derive = (S00abs - thrust::abs(S00(sm-EPS)) )/EPS;
    //fptype s00sabs = thrust::abs(S00(s));
    
    fptype func=0.;

    if(ni==1){
        func= exp(- pow( (sqrt(-log(S00abs)) + (-8.*q(sm)*(q(s) - q(sm))*deri2)/(2.*sqrt(-log(S00abs))*S00abs) + e00*pow((q(s)/q(sm) - 1.),2.) + e01*pow((q(s)/q(sm) - 1.),3.)),2.));
        //printf("S00abs = %f \n",S00abs);
    }
    if(ni==2){
        func=exp(- pow( (sqrt(-log(S00abs)) + (-8.*q(sm)*(q(s) - q(sm))*deri2p)/(2.*sqrt(-log(S00abs))*S00abs) + e00*pow((q(s)/q(sm) - 1.),2.) + e01*pow((q(s)/q(sm) - 1.),3.) + e02*pow((q(s)/q(sm) - 1.),4.)),2.));
    }
    if(ni==3){
        func= exp(- pow( (sqrt(-log(S00abs)) + (-8.*q(sm)*(q(s) - q(sm))*deri2pp)/(2.*sqrt(-log(S00abs))*S00abs)  + e00*pow((q(s)/q(sm) - 1.),2.) + e01*pow((q(s)/q(sm) - 1.),3.) + e02*pow((q(s)/q(sm) - 1.),4.)),2.));
    }
    if(s<sm){ // s<1.96GeV^2
        S00abs = thrust::abs(S00(s));
        return S00abs;
    } else return func;
}

__device__ fpcomplex ampt00 (fptype s){
    fpcomplex eta(Inela(s),0.);
    fpcomplex c(0.,2.*sqrt( 1.0 - (4.*mp*mp)/s));
    fpcomplex argpp(0.,2.*argument(s)*M_PI/180.0);
    return (eta*thrust::exp(argpp) - fpcomplex(1.,0.))/c ;
}


} // namespace GooFit
