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
        q = 1.0;
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
            p = 1.0;
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
            pstar = 1.0;
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
    fptype m12,
    fptype m13,
    fptype m23,
    unsigned int cyclic_index) -> fptype {

        auto const _mA  = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
        auto const _mC  = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
        auto const _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12));
        auto const _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    

        fptype EACmsAB = (_mAB - _mA*_mA + _mC*_mC)/(2.0*sqrt(_mAB));
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

        if(cyclic_index==PAIR_12 || cyclic_index==PAIR_13)
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

} // namespace GooFit
