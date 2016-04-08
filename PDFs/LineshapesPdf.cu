/*
04/05/2016
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

TODO:
- reorganize this file into multiple files so lineshapes and spinfactors are seperated. Also all helper functions should go into seperate file.

This file includes some lineshapes and spinfactors.
Also right now it is the home to some helper functions needed and an implementation of a simple 4-vec class that works on the GPU
*/



#include "LineshapesPdf.hh" 



EXEC_TARGET fptype Mass(const fptype* P0){
  return SQRT(-P0[0]*P0[0] - P0[1]*P0[1] - P0[2]*P0[2] + P0[3]*P0[3]);
}
EXEC_TARGET fptype Mass(const fptype* P0, const fptype* P1){
  return SQRT( -((P0[0]+P1[0]) * (P0[0]+P1[0])) - ((P0[1]+P1[1]) * (P0[1]+P1[1])) - ((P0[2]+P1[2]) * (P0[2]+P1[2])) + ((P0[3]+P1[3]) * (P0[3]+P1[3])) );
}
EXEC_TARGET fptype Mass(const fptype* P0, const fptype* P1, const fptype* P2){
  return SQRT( -((P0[0]+P1[0]+P2[0]) * (P0[0]+P1[0]+P2[0])) - ((P0[1]+P1[1]+P2[1]) * (P0[1]+P1[1]+P2[1])) - ((P0[2]+P1[2]+P2[2]) * (P0[2]+P1[2]+P2[2])) + ((P0[3]+P1[3]+P2[3]) * (P0[3]+P1[3]+P2[3])) );
}
EXEC_TARGET fptype VecDot(const fptype* P0, const fptype* P1){
  return ( P0[0]*P1[0]  + P0[1]+P1[1] + P0[2]+P1[2] + P0[3]+P1[3] );
}


EXEC_TARGET gpuLVec::gpuLVec(fptype x, fptype y, fptype z, fptype e)  :  X(x), Y(y), Z(z), E(e){};
EXEC_TARGET fptype gpuLVec::Dot(const gpuLVec& rhs) const { return E*rhs.E - X*rhs.X - Y*rhs.Y - Z*rhs.Z ;}
EXEC_TARGET gpuLVec& gpuLVec::operator+=(const gpuLVec& rhs){
  X+=rhs.X;
  Y+=rhs.Y;
  Z+=rhs.Z;
  E+=rhs.E;
  return *this;
}                          
EXEC_TARGET gpuLVec& gpuLVec::operator-=(const gpuLVec& rhs){
  X-=rhs.X;
  Y-=rhs.Y;
  Z-=rhs.Z;
  E-=rhs.E;
  return *this;
}    
EXEC_TARGET gpuLVec operator+(gpuLVec lhs, const gpuLVec& rhs){return lhs+=rhs;}
EXEC_TARGET gpuLVec operator-(gpuLVec lhs, const gpuLVec& rhs){return lhs-=rhs;}


EXEC_TARGET fptype LeviCevita(const gpuLVec& p1, const gpuLVec& p2, const gpuLVec& p3, const gpuLVec& p4){
  // this calculates the determinant of the 4x4 matrix build out of p1,p2,p3,p4
  return
     p1.getZ() * p2.getY() * p3.getX() * p4.getE() - p1.getY() * p2.getZ() * p3.getX() * p4.getE() -
     p1.getZ() * p2.getX() * p3.getY() * p4.getE() + p1.getX() * p2.getZ() * p3.getY() * p4.getE() +
     p1.getY() * p2.getX() * p3.getZ() * p4.getE() - p1.getX() * p2.getY() * p3.getZ() * p4.getE() -
     p1.getZ() * p2.getY() * p3.getE() * p4.getX() + p1.getY() * p2.getZ() * p3.getE() * p4.getX() +
     p1.getZ() * p2.getE() * p3.getY() * p4.getX() - p1.getE() * p2.getZ() * p3.getY() * p4.getX() -
     p1.getY() * p2.getE() * p3.getZ() * p4.getX() + p1.getE() * p2.getY() * p3.getZ() * p4.getX() +
     p1.getZ() * p2.getX() * p3.getE() * p4.getY() - p1.getX() * p2.getZ() * p3.getE() * p4.getY() -
     p1.getZ() * p2.getE() * p3.getX() * p4.getY() + p1.getE() * p2.getZ() * p3.getX() * p4.getY() +
     p1.getX() * p2.getE() * p3.getZ() * p4.getY() - p1.getE() * p2.getX() * p3.getZ() * p4.getY() -
     p1.getY() * p2.getX() * p3.getE() * p4.getZ() + p1.getX() * p2.getY() * p3.getE() * p4.getZ() +
     p1.getY() * p2.getE() * p3.getX() * p4.getZ() - p1.getE() * p2.getY() * p3.getX() * p4.getZ() -
     p1.getX() * p2.getE() * p3.getY() * p4.getZ() + p1.getE() * p2.getX() * p3.getY() * p4.getZ();
}


EXEC_TARGET fptype S_VV_PPPP_S (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",0, P1.getX(), P1.getY(), P1.getZ(),P1.getE());
  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",1, P2.getX(), P2.getY(), P2.getZ(),P2.getE());
  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",2, P3.getX(), P3.getY(), P3.getZ(),P3.getE());
  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",3, P4.getX(), P4.getY(), P4.getZ(),P4.getE());

  gpuLVec pV1 = P1 + P2;
  gpuLVec qV1 = P1 - P2;
  gpuLVec pV2 = P3 + P4;
  gpuLVec qV2 = P3 - P4;
  
  fptype MV1 = SQRT(pV1.Dot(pV1));
  fptype MV2 = SQRT(pV2.Dot(pV2));

  fptype returnVal = (qV1.Dot(qV2) 
                   - qV1.Dot(pV1) * pV1.Dot(qV2) / (MV1*MV1)
                   - qV1.Dot(pV2) * pV2.Dot(qV2) / (MV2*MV2)
                   + qV1.Dot(pV1) * pV1.Dot(pV2) * pV2.Dot(qV2) 
                   / (MV1*MV1 * MV2*MV2));
  // printf("s1 %.5g; %i,%i,%i,%i\n",returnVal, indices[2], indices[3], indices[4], indices[5]);
  return returnVal;
}

EXEC_TARGET fptype S_VV_PPPP_P (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV1 = P1 + P2;
  gpuLVec qV1 = P1 - P2;
  gpuLVec pV2 = P3 + P4;
  gpuLVec qV2 = P3 - P4;
  
  gpuLVec pD = pV1 + pV2;
  gpuLVec qD = pV1 - pV2;

  return LeviCevita(pD, qD, qV1, qV2);
}

EXEC_TARGET fptype S_VV_PPPP_D (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV1 = P1 + P2;
  gpuLVec qV1 = P1 - P2;
  gpuLVec pV2 = P3 + P4;
  gpuLVec qV2 = P3 - P4;
  
  // printf("%f, %f, %f, %f\n",P1.getX(), P1.getY(), P1.getZ(), P1.getE() );
  // printf("%f, %f, %f, %f\n",P2.getX(), P2.getY(), P2.getZ(), P2.getE() );
  // printf("%f, %f, %f, %f\n",P3.getX(), P3.getY(), P3.getZ(), P3.getE() );
  // printf("%f, %f, %f, %f\n",P4.getX(), P4.getY(), P4.getZ(), P4.getE() );

  fptype MV1 = SQRT(pV1.Dot(pV1));
  fptype MV2 = SQRT(pV2.Dot(pV2));
  fptype returnVal = (  qV1.Dot(pV2) - qV1.Dot(pV1) * pV1.Dot(pV2)/(MV1*MV1)
                     )*( 
                     qV2.Dot(pV1) - qV2.Dot(pV2) * pV2.Dot(pV1)/(MV2*MV2)
                     );
  return returnVal;
}

EXEC_TARGET fptype S_AP1_AtoVP2_VtoP3P4 (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV = P3 + P4;
  gpuLVec qV = P3 - P4;
  gpuLVec pA = P2 + pV;
  gpuLVec p0 = P1;  
  gpuLVec pD = P1 + pA;
  gpuLVec qD = pA - P1;
  
  fptype MA = SQRT(pA.Dot(pA));
  fptype MV = SQRT(pV.Dot(pV));
  fptype returnVal =  P1.Dot(qV)
      -   p0.Dot(pA) * pA.Dot(qV) / (MA*MA)
      -   p0.Dot(pV) * pV.Dot(qV) / (MV*MV)
      +   p0.Dot(pA) * pA.Dot(pV) * pV.Dot(qV) / (MA*MA * MV*MV);
  // printf("spin %.7g\n",returnVal );
  return returnVal;
}



EXEC_TARGET devcomplex<fptype> BW_DP (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+0];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital          = indices[4];

  fptype frFactor = 1;

  fptype rMass2 = Mpair * Mpair;
  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame 
  fptype measureDaughterMoms = twoBodyCMmom(rMass2, m1, m2);
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, m1, m2);

  if (0 != orbital) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, orbital, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, orbital, meson_radius); 
  }  
  // RBW evaluation
  fptype A = (resmass - rMass2); 
  fptype B = resmass*reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*orbital + 1) * frFactor / SQRT(rMass2);
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> ret(A*C, B*C); // Dropping F_D=1

  ret *= SQRT(frFactor); 
  return ret; 
}

EXEC_TARGET devcomplex<fptype> BW_MINT (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+0];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital          = indices[4];


  const unsigned int to2Lplus1    = 2 * orbital + 1;

  fptype mass = resmass;
  fptype width = reswidth;
  fptype mumsRecoMass2 = Mpair*Mpair;
  

  fptype mpsq = (m1+m2)*(m1+m2);
  fptype mmsq = (m1-m2)*(m1-m2);
  fptype num  = (mumsRecoMass2 - mpsq)*(mumsRecoMass2 - mmsq);
  fptype num2  = (mass*mass - mpsq)*(mass*mass - mmsq);
  fptype pABSq = num/(4*mumsRecoMass2);
  fptype prSqForGofM = num2/(4*mass*mass);
  fptype pratio = SQRT(pABSq/prSqForGofM);

  fptype pratio_to_2Jplus1 = 1;

  for(int i=0; i < to2Lplus1; i++){
    pratio_to_2Jplus1 *= pratio;
  }    

  fptype mratio = mass/Mpair;

  fptype thisFR = SQRT((1+meson_radius*meson_radius*prSqForGofM) / (1+meson_radius*meson_radius*pABSq));

  fptype GofM = width * pratio_to_2Jplus1 *mratio * thisFR * thisFR;

  fptype gamma = SQRT(mass*mass*(mass*mass + width*width));
  fptype k     = mass*width*gamma/SQRT(mass*mass+gamma);

  devcomplex<fptype> BW(mass*mass - mumsRecoMass2, mass*GofM);
  fptype den = (mass*mass - mumsRecoMass2) * (mass*mass - mumsRecoMass2) + mass * GofM * mass * GofM;

  devcomplex<fptype> ret = (SQRT(k) * thisFR)/den * BW;

  // printf("%.7g, %.7g, %.7g, %i, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g\n", m1, m2, Mpair, to2Lplus1, GofM, pratio_to_2Jplus1, mratio, k , meson_radius, prSqForGofM, thisFR, ret.real, ret.imag );
  return  ret ; 
}

EXEC_TARGET devcomplex<fptype> bugg_rho2(const fptype& s, const fptype m){
  fptype rho_squared = 1. - 4. * m*m /s;
  devcomplex<fptype> returnVal = (rho_squared >= 0) ? devcomplex<fptype>(1,0) : devcomplex<fptype>(0,1);
  rho_squared = (rho_squared >= 0) ? SQRT(rho_squared) : SQRT(-rho_squared);
  return rho_squared * returnVal;
}

EXEC_TARGET fptype bugg_j1(const fptype& s, const fptype m){
  fptype rho_pipi = bugg_rho2(s, m).real;
  fptype returnVal = 2.;
  returnVal += (rho_pipi>0.) ? rho_pipi * LOG((1.-rho_pipi)/(1.+rho_pipi)) : 0;
  return returnVal/M_PI;
}

EXEC_TARGET fptype bugg_Gamma_4pi(const fptype& s, const fptype mpi, const fptype& g_4pi, const fptype& M, const fptype& lambda_4pi, const fptype& s0_4pi){
  fptype returnVal = (s < (16. * mpi*mpi)) ? 0 : g_4pi* (1./(1+EXP(lambda_4pi*(s0_4pi-s))))/(1./(1+EXP(lambda_4pi*(s0_4pi-M*M))));
  return returnVal;
}

// this lineshape is not tested yet!
EXEC_TARGET devcomplex<fptype> bugg_MINT (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype resmass                = cudaArray[indices[2]];
  // fptype reswidth               = cudaArray[indices[3]];
  unsigned int spin             = indices[4];
  fptype frFactor               = 1;
  fptype s                      = Mpair*Mpair;
  resmass                      *= resmass;
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2).
  
  fptype measureDaughterMoms = twoBodyCMmom(s, m1, m2);
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, m1, m2);

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius);
  }
  fptype M            = 0.935;
  fptype b1           = 1.302;
  fptype b2           = 0.340;
  fptype A            = 2.426;
  fptype g_4pi        = 0.011;
  fptype g_2K         = 0.6;
  fptype g_2eta       = 0.2;
  fptype alpha        = 1.3;
  fptype sA           = 0.41;
  fptype s0_4pi       = 7.082/2.845;
  fptype lambda_4pi   = 2.845;
  fptype mPiPlus      = .13957018;
  fptype mKPlus       = .493677;
  fptype mEta         = .54751;

  fptype g1sq = (b1+b2*s)*EXP(-(s-M*M)/A);
  fptype z = bugg_j1(s, mPiPlus) - bugg_j1(M*M, mPiPlus);


  devcomplex<fptype> gamma_2pi = devcomplex<fptype>(g1sq * (s -sA*mPiPlus*mPiPlus)/(M*M -sA*mPiPlus*mPiPlus)*bugg_rho2(s, mPiPlus).real, 0);
  devcomplex<fptype> gamma_2K = g_2K * g1sq * s/(M*M) * EXP((-1)*alpha * SQRT((s-4.*mKPlus*mKPlus)*(s-4.*mKPlus*mKPlus))) * bugg_rho2(s, mKPlus);
  devcomplex<fptype> gamma_2eta = g_2eta * g1sq * s/(M*M) * EXP((-1)*alpha * SQRT((s-4.*mEta*mEta)*(s-4.*mEta*mEta))) * bugg_rho2(s, mEta);
  devcomplex<fptype> gamma_4pi = devcomplex<fptype>(bugg_Gamma_4pi(s,mPiPlus,g_4pi,M,lambda_4pi,s0_4pi),0);
  
  devcomplex<fptype> Gamma_tot = gamma_2pi + gamma_2K + gamma_2eta + gamma_4pi;
  
  devcomplex<fptype> num = M * gamma_2pi;
  devcomplex<fptype> den = devcomplex<fptype>(M*M - s - g1sq * (s-sA*mPiPlus*mPiPlus) / (M*M-sA*mPiPlus*mPiPlus) * z,0) - devcomplex<fptype>(0,1) * M * Gamma_tot;

  return num/den;
}



// this lineshape is not tested yet!
EXEC_TARGET devcomplex<fptype> lass_DP (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital             = indices[4];
  fptype frFactor               = 1;
  fptype rMass2                 = Mpair*Mpair;
  resmass                      *= resmass;
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2).
  
  fptype measureDaughterMoms = twoBodyCMmom(rMass2, m1, m2);
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, m1, m2);

  if (0 != orbital) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, orbital, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, orbital, meson_radius);
  }

  //Implement LASS:
  /*
  fptype s = kinematics(m12, m13, _trackinfo[i]);
  fptype q = twoBodyCMmom(s, _trackinfo[i]);
  fptype m0  = _massRes[i]->getValFast();
  fptype _g0 = _gammaRes[i]->getValFast();
  int spin   = _spinRes[i];
  fptype g = runningWidthFast(s, m0, _g0, spin, _trackinfo[i], FrEval(s, m0, _trackinfo[i], spin));
  */

  fptype q = measureDaughterMoms;
  fptype g = reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*orbital + 1) * frFactor / SQRT(rMass2);

  fptype _a    = 0.22357;
  fptype _r    = -15.042;
  fptype _R    = 1; // ?
  fptype _phiR = 1.10644;
  fptype _B    = 0.614463;
  fptype _phiB = -0.0981907;

  // background phase motion
  fptype cot_deltaB = (1.0 / (_a*q)) + 0.5*_r*q;
  fptype qcot_deltaB = (1.0 / _a) + 0.5*_r*q*q;

  // calculate resonant part
  devcomplex<fptype> expi2deltaB = devcomplex<fptype>(qcot_deltaB,q)/devcomplex<fptype>(qcot_deltaB,-q);
  devcomplex<fptype>  resT = devcomplex<fptype>(cos(_phiR+2*_phiB),sin(_phiR+2*_phiB))*_R;

  devcomplex<fptype> prop = devcomplex<fptype>(1, 0)/devcomplex<fptype>(resmass-rMass2, SQRT(resmass)*g);
  // resT *= prop*m0*_g0*m0/twoBodyCMmom(m0*m0, _trackinfo[i])*expi2deltaB;
  resT *= prop*(resmass*reswidth/nominalDaughterMoms)*expi2deltaB;

  // calculate bkg part
  resT += devcomplex<fptype>(cos(_phiB),sin(_phiB))*_B*(cos(_phiB)+cot_deltaB*sin(_phiB))*SQRT(rMass2)/devcomplex<fptype>(qcot_deltaB,-q);

  resT *= SQRT(frFactor);
  return resT;
}


EXEC_TARGET devcomplex<fptype> nonres_DP (fptype m12, fptype m1, fptype m2, unsigned int* indices) {
  return devcomplex<fptype>(1, 0); 
}

MEM_DEVICE resonance_function_ptr ptr_to_BW_DP = BW_DP;
MEM_DEVICE resonance_function_ptr ptr_to_BW_MINT = BW_MINT;
MEM_DEVICE resonance_function_ptr ptr_to_lass_DP = lass_DP;
MEM_DEVICE resonance_function_ptr ptr_to_NONRES_DP = nonres_DP;
MEM_DEVICE spin_function_ptr ptr_to_S_VV_PPPP_S = S_VV_PPPP_S;
MEM_DEVICE spin_function_ptr ptr_to_S_VV_PPPP_P = S_VV_PPPP_P;
MEM_DEVICE spin_function_ptr ptr_to_S_VV_PPPP_D = S_VV_PPPP_D;
MEM_DEVICE spin_function_ptr ptr_to_S_AP1_AtoVP2_VtoP3P4 = S_AP1_AtoVP2_VtoP3P4
;


Lineshape::Lineshape (string name,
						Variable* mass, 
						Variable* width, 
						unsigned int L, 
						unsigned int cyc,
            bool useMINTBW = false) 
  : GooPdf(0, name)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Making room for index of decay-related constants. Assumption:
  // These are mother mass and three daughter masses in that order.
  // They will be registered by the object that uses this resonance,
  // which will tell this object where to find them by calling setConstantIndex. 
  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(L);
  pindices.push_back(cyc); 

  if(useMINTBW) GET_FUNCTION_ADDR(ptr_to_BW_MINT);
  else GET_FUNCTION_ADDR(ptr_to_BW_DP);
  initialise(pindices); 
}



Lineshape::Lineshape (string name)
  : GooPdf(0, name)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  GET_FUNCTION_ADDR(ptr_to_NONRES_DP);
  initialise(pindices); 
}

Amplitude::Amplitude(std::string uniqueDecayStr, Variable* ar, Variable* ai, std::map<std::string, Lineshape*> LS, std::map<std::string, SpinFactor*> SF, unsigned int nPerm)
  : _uniqueDecayStr(uniqueDecayStr),
    _ar(ar),
    _ai(ai),
    _LS(LS),
    _SF(SF),
    _nPerm(nPerm)
{}



SpinFactor::SpinFactor (std::string name, unsigned int kind, unsigned int P0, unsigned int P1, unsigned int P2, unsigned int P3)
 : GooPdf(0,name){
  vector<unsigned int> pindices; 
  pindices.push_back(0); //dummy for index to constants.
  pindices.push_back(P0);
  pindices.push_back(P1);
  pindices.push_back(P2);
  pindices.push_back(P3);
  switch(kind){
    case 0:
      GET_FUNCTION_ADDR(ptr_to_S_VV_PPPP_S);
      break;
    case 1:
      GET_FUNCTION_ADDR(ptr_to_S_VV_PPPP_P);
      break;
    case 2:
      GET_FUNCTION_ADDR(ptr_to_S_VV_PPPP_D);
      break;
    case 3:
      GET_FUNCTION_ADDR(ptr_to_S_AP1_AtoVP2_VtoP3P4);
      break;
    
    default:
      std::cout << "No Spinfunction implemented for that kind." << std::endl;
      exit(0);
      break;
  }
  
  initialise(pindices);
}


EXEC_TARGET void get4Vecs (fptype* Vecs, const unsigned int& constants, const fptype& m12, const fptype& m34, const fptype& cos12, const fptype& cos34, const fptype& phi){
  fptype M = functorConstants[constants + 1]; 
  fptype m1  = functorConstants[constants + 2]; 
  fptype m2  = functorConstants[constants + 3]; 
  fptype m3  = functorConstants[constants + 4]; 
  fptype m4  = functorConstants[constants + 5]; 
  // printf("g4v %f, %f, %f, %f, %f\n",M, m1, m2, m3, m4 );
  fptype E1 = (m12*m12 + m1*m1 - m2*m2) / (2 * m12) ; 
  fptype E2 = (m12*m12 - m1*m1 + m2*m2) / (2 * m12) ; 
  fptype E3 = (m34*m34 + m3*m3 - m4*m4) / (2 * m34) ; 
  fptype E4 = (m34*m34 - m3*m3 + m4*m4) / (2 * m34) ; 
  fptype p1 = SQRT(E1*E1 - m1*m1);
  fptype p3 = SQRT(E3*E3 - m3*m3); 
  fptype sin12 = SQRT(1-cos12*cos12);
  fptype sin34 = SQRT(1-cos34*cos34);
  fptype ED1 = ( M*M + m12*m12 - m34*m34) / (2*m12);
  fptype PD1 = SQRT(ED1*ED1 - M*M);
  fptype beta1 = PD1 / ED1;
  fptype gamma1 = 1.0/SQRT(1-beta1*beta1);
  fptype ED2 = ( M*M - m12*m12 + m34*m34) / (2*m34);
  fptype PD2 = SQRT(ED2*ED2 - M*M);
  fptype beta2 = -PD2 / ED2;
  fptype gamma2 = 1.0/SQRT(1-beta2*beta2);
  // printf("g4v %f, %f, %f, %f, %f\n",E1, m1, E2, p1, p3 );

  //set X-component
  Vecs[0] = cos12*p1;
  Vecs[4] = -cos12*p1;
  Vecs[8] = -cos34*p3;
  Vecs[12] = cos34*p3;
 
  //set Y-component
  Vecs[1] = sin12*p1;
  Vecs[5] = -sin12*p1;
  Vecs[9] = -sin34*p3;
  Vecs[13] = sin34*p3;

  //set Y-component
  Vecs[2]  = 0;
  Vecs[6]  = 0;
  Vecs[10] = 0;
  Vecs[14] = 0;

  //set E-component
  Vecs[3]  = E1;
  Vecs[7]  = E2;
  Vecs[11] = E3;
  Vecs[15] = E4;

  fptype tmpE = Vecs[3];
  fptype tmpX = Vecs[0];
  Vecs[3] = gamma1*( tmpE + beta1*tmpX );
  Vecs[0] = gamma1*( tmpX + beta1*tmpE );

  tmpE = Vecs[7];
   tmpX = Vecs[4];
  Vecs[7] = gamma1*( tmpE + beta1*tmpX );
  Vecs[4] = gamma1*( tmpX + beta1*tmpE );

  tmpE = Vecs[11];
  tmpX = Vecs[8];
  Vecs[11] = gamma2*( tmpE + beta2*tmpX );
  Vecs[8] = gamma2*( tmpX + beta2*tmpE );

  tmpE = Vecs[15];
  tmpX = Vecs[12];
  Vecs[15] = gamma2*( tmpE + beta2*tmpX );
  Vecs[12] = gamma2*( tmpX + beta2*tmpE );

  // rotation around X-axis of the first two vectors.
  fptype cosphi = cos(phi);
  fptype sinphi = sin(phi);

  // note that Z-component is zero thus rotation is as easy as this:
  Vecs[2] = sinphi*Vecs[1]; 
  Vecs[1] = cosphi*Vecs[1];

  Vecs[6] = sinphi*Vecs[5]; 
  Vecs[5] = cosphi*Vecs[5];
  
} 

EXEC_TARGET fptype getmass(const unsigned int& pair, fptype& d1, fptype& d2, const fptype* vecs, const fptype& m1, const fptype& m2, const fptype& m3, const fptype& m4){
    const fptype* P1 = vecs;
    const fptype* P2 = (vecs+4);
    const fptype* P3 = (vecs+8);
    const fptype* P4 = (vecs+12);
    fptype mpair;
  switch(pair){

    case 2:
      d1 = m1;
      d2 = m3;
      mpair = Mass(P1,P3);
    break;

    case 3:
      d1 = m1;
      d2 = m4;
      mpair = Mass(P1,P4);
    break;

    case 4:
      d1 = m2;
      d2 = m3;
      mpair = Mass(P3,P3);
    break;
    
    case 5:
      d1 = m2;
      d2 = m4;
      mpair = Mass(P2,P4);
    break;

    case 6:
      d1 = Mass(P1,P2);
      d2 = m3;
      mpair = Mass(P1,P2,P3);
    break;

    case 7:
      d1 = Mass(P1,P3);
      d2 = m2;
      mpair = Mass(P1,P2,P3);
    break;

    case 8:
      d1 = Mass(P2,P3);
      d2 = m1;
      mpair = Mass(P1,P2,P3);
    break;

    case 9:
      d1 = Mass(P1,P2);
      d2 = m4;
      mpair = Mass(P1,P2,P4);
    break;

    case 10:
      d1 = Mass(P1,P4);
      d2 = m2;
      mpair = Mass(P1,P2,P4);
    break;

    case 11:
      d1 = Mass(P2,P4);
      d2 = m1;
      mpair = Mass(P1,P2,P4);
    break;

    case 12:
      d1 = Mass(P1,P3);
      d2 = m4;
      mpair = Mass(P1,P3,P4);
    break;

    case 13:
      d1 = Mass(P1,P4);
      d2 = m3;
      mpair = Mass(P1,P3,P4);
    break;

    case 14:
      d1 = Mass(P3,P4);
      d2 = m1;
      mpair = Mass(P1,P3,P4);
    break;

    case 15:
      d1 = Mass(P2,P3);
      d2 = m4;
      mpair = Mass(P2,P3,P4);
    break;

    case 16:
      d1 = Mass(P2,P4);
      d2 = m3;
      mpair = Mass(P2,P3,P4);
    break;

    case 17:
      d1 = Mass(P3,P4);
      d2 = m2;
      mpair = Mass(P2,P3,P4);
    break;
  }
  return mpair;
}

