/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

This file includes some lineshapes and spinfactors.
Also right now it is the home to some helper functions needed and an implementation of a simple 4-vec class that works on the GPU
*/

#include "LineshapesPdf.hh" 


// Form factors as in pdg http://pdg.lbl.gov/2012/reviews/rpp2012-rev-dalitz-analysis-formalism.pdf
EXEC_TARGET fptype BL_PRIME (fptype z2, fptype z02, int L) {
  if(0 ==L) return 1.0;
  else if (1 == L) return(1+z02)/(1+z2);
  else if (2 == L) return ( z02*z02 + 3*z02 + 9 ) / ( z2*z2 + 3*z2 + 9 )  ;
  else{
    printf("ERROR! Oribtal > 2 not supported!\n");
    return 0;
  }
  // Spin 3 and up not accounted for. 
}

EXEC_TARGET fptype BL(fptype z2, int L) {
  if(0 ==L) return 1.0;
  else if( 1==L) return 2*z2/(1+z2);
  else if( 2==L) return (13*z2*z2 )/ ( z2*z2 + 3*z2 + 9 )  ;
  else{
    printf("ERROR! Oribtal > 2 not supported!\n");
    return 0;
  }
  // Spin 3 and up not accounted for. 
}

//This function is modeled after BW_BW::getVal() in BW_BW.cpp from the MINT package written by Jonas Rademacker. 
EXEC_TARGET devcomplex<fptype> BW (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+0];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital          = indices[4];
  unsigned int FF               = indices[6];

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
  fptype prSq2 = prSqForGofM < 0 ? 0 : prSqForGofM;
  prSqForGofM = FABS(prSqForGofM);

  fptype pratio = SQRT(pABSq/prSqForGofM);

  fptype pratio_to_2Jplus1 = 1;

  for(int i=0; i < to2Lplus1; i++){
    pratio_to_2Jplus1 *= pratio;
  }    

  fptype mratio = mass/Mpair;
  fptype r2 = meson_radius*meson_radius;
  fptype thisFR = BL_PRIME(pABSq*r2, prSqForGofM*r2, orbital);
  fptype frFactor = 1;
  if (0 != orbital and 0 != FF) {
    frFactor =  (FF==1? BL(pABSq*r2, orbital) : BL_PRIME(pABSq*r2, prSq2*r2, orbital));
  } 

  fptype GofM = width * pratio_to_2Jplus1 *mratio * thisFR;

  fptype gamma = SQRT(mass*mass*(mass*mass + width*width));
  fptype k     = mass*width*gamma/SQRT(mass*mass+gamma);

  devcomplex<fptype> BW(mass*mass - mumsRecoMass2, mass*GofM);
  fptype den = (mass*mass - mumsRecoMass2) * (mass*mass - mumsRecoMass2) + mass * GofM * mass * GofM;

  devcomplex<fptype> ret = (SQRT(k * frFactor))/den * BW;
  // printf("m1, m2, Mpair, to2Lplus1, GofM, thisFR, pratio, mratio, pABSq , prSqForGofM, FF, ret.real, ret.imag\n");
  // printf("BW %.7g, %.7g, %.7g, %i, %i, %i, %i\n",meson_radius, resmass, reswidth, orbital, FF, indices[2], indices[3]);
  // printf("BW %.7g, %.7g, %.7g, %i, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g\n", m1, m2, Mpair, to2Lplus1, GofM, thisFR, pratio, mratio, pABSq, prSqForGofM, frFactor, ret.real, ret.imag );
  return  ret ; 
}

//This function is modeled after SBW from the MINT package written by Jonas Rademacker. 
EXEC_TARGET devcomplex<fptype> SBW (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+0];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital          = indices[4];
  unsigned int FF               = indices[6];

  fptype mass = resmass;
  fptype width = reswidth;
  fptype mumsRecoMass2 = Mpair*Mpair;
  

  fptype mpsq = (m1+m2)*(m1+m2);
  fptype mmsq = (m1-m2)*(m1-m2);
  fptype num  = (mumsRecoMass2 - mpsq)*(mumsRecoMass2 - mmsq);
  fptype num2  = (mass*mass - mpsq)*(mass*mass - mmsq);
  fptype pABSq = num/(4*mumsRecoMass2);
  fptype prSq = num2/(4*mass*mass);
  fptype prSq2 = prSq < 0 ? 0 : prSq;
  prSq = FABS(prSq);

  fptype r2 = meson_radius*meson_radius;
  fptype frFactor = 1;
  if (0 != orbital and 0 != FF) {
    frFactor =  (FF==1? BL(pABSq*r2, orbital) : BL_PRIME(pABSq*r2, prSq2*r2, orbital));
  } 

  fptype GofM = width;

  fptype gamma = SQRT(mass*mass*(mass*mass + width*width));
  fptype k     = mass*width*gamma/SQRT(mass*mass+gamma);

  devcomplex<fptype> BW(mass*mass - mumsRecoMass2, mass*GofM);
  fptype den = (mass*mass - mumsRecoMass2) * (mass*mass - mumsRecoMass2) + mass * GofM * mass * GofM;

  devcomplex<fptype> ret = (SQRT(k * frFactor))/den * BW;

  // printf("m1, m2, Mpair, GofM, pABSq , prSq, FF, ret.real, ret.imag\n");
  // printf("SBW %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g\n", m1, m2, Mpair, GofM, pABSq, prSq, frFactor, ret.real, ret.imag );
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

//This function is an adaptation from the bugg lineshape implemented in the MINT package written by Jonas Rademacker. 
// this lineshape is not tested yet!
EXEC_TARGET devcomplex<fptype> bugg_MINT (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype s                      = Mpair*Mpair;

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
  
  // devcomplex<fptype> num = M * gamma_2pi; //only for elastic scattering, not production 
  devcomplex<fptype> den = devcomplex<fptype>(M*M - s - M * g1sq * (s-sA*mPiPlus*mPiPlus) / (M*M-sA*mPiPlus*mPiPlus) * z,0) - devcomplex<fptype>(0,1) * M * Gamma_tot;
  devcomplex<fptype> returnVal = 1.0/den;
  // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",gamma_2pi.real, gamma_2pi.imag, gamma_2K.real, gamma_2K.imag, gamma_2eta.real, gamma_2eta.imag, gamma_4pi.real, gamma_4pi.imag);
  // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Mpair, Gamma_tot.real, Gamma_tot.imag, g1sq, z, den.real, den.imag, returnVal.real, returnVal.imag);
  
  //the factor sqrt(1000) gives the correct result in comparison with mint2, I think its because BW/SBW 
  // have a factor of sqrt(k) which these lineshapes dont have. For now this stays because it works. further investigation needed.
  return returnVal*SQRT(1000.0);
}

EXEC_TARGET devcomplex<fptype> lass_MINT (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital             = indices[4];
  fptype frFactor               = 1;
  fptype rMass2                 = Mpair*Mpair;
  
  fptype a = 2.07;
  fptype r = 3.32;
  fptype phi = 0.0;
  fptype cutCoff = 1.8;

  fptype mpsq = (m1+m2)*(m1+m2);
  fptype mmsq = (m1-m2)*(m1-m2);
  fptype num  = (rMass2 - mpsq)*(rMass2 - mmsq);
  fptype num2  = (resmass*resmass - mpsq)*(resmass*resmass - mmsq);
  fptype pABSq = num/(4*rMass2);
  fptype prSq = FABS(num2/(4*resmass*resmass));

  fptype y = 2.0 * a*SQRT(pABSq);
  fptype x = 2.0 + a * r * pABSq;
  fptype cotDeltaBg = x / y;
  devcomplex<fptype> phaseshift((cotDeltaBg*cotDeltaBg-1)/(1+cotDeltaBg*cotDeltaBg), 2*cotDeltaBg / ( 1 + cotDeltaBg*cotDeltaBg));
  // (cotDeltaBg*cotDeltaBg-1)/(1+cotDeltaBg*cotDeltaBg) = cos(2*delta)     2*cotDeltaBg / ( 1 + cotDeltaBg*cotDeltaBg) = sin(2*delta)
  devcomplex<fptype> den(SQRT(pABSq)*cotDeltaBg, (-1.)*SQRT(pABSq));
  fptype SF = Mpair * SQRT(prSq) / ( resmass * resmass*reswidth );
  devcomplex<fptype> BG = SF / den ;
  devcomplex<fptype> returnVal = BG + phaseshift * BW(Mpair, m1, m2, indices);
  // printf("Lass: %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",phaseshift.real, phaseshift.imag, den.real, den.imag, BG.real, BG.imag,  returnVal.real, returnVal.imag);
  
  return returnVal;
}

//Lass lineshape as implemented in MINT3 by Tim Evans. Not quite sure why this is different than the Mint2 version.
EXEC_TARGET devcomplex<fptype> lass_MINT3 (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital             = indices[4];
  fptype frFactor               = 1;
  fptype rMass2                 = Mpair*Mpair;
  
  fptype a = 2.07;
  fptype r = 3.32;

  fptype mpsq = (m1+m2)*(m1+m2);
  fptype mmsq = (m1-m2)*(m1-m2);
  fptype num  = (rMass2 - mpsq)*(rMass2 - mmsq);
  fptype num2  = (resmass*resmass - mpsq)*(resmass*resmass - mmsq);
  fptype pABSq = num/(4*rMass2);
  fptype prSq = FABS(num2/(4*resmass*resmass));

  fptype pratio = SQRT(pABSq/prSq);

  fptype pratio_to_2Jplus1 = 1;

  for(int i=0; i < 2*orbital+1; i++){
    pratio_to_2Jplus1 *= pratio;
  }

  fptype mratio = resmass/Mpair;
  fptype r2 = meson_radius*meson_radius;
  fptype thisFR = BL_PRIME(pABSq*r2, prSq*r2, orbital);
  fptype GofM = reswidth * pratio_to_2Jplus1 *mratio * thisFR;


  fptype y = 2.0 * a*SQRT(pABSq);
  fptype x = 2.0 + a * r * pABSq;
  fptype phase = atan(y/x) +  atan(resmass*GofM/(resmass*resmass - rMass2)) ;
  fptype rho = 1.0 / SQRT(pABSq/rMass2);
  devcomplex<fptype> returnVal = SIN(phase) * devcomplex<fptype>(COS(phase),SIN(phase)) * rho;
  // printf("Lass: %.5g %.5g %.5g\n",returnVal.real, returnVal.imag, GofM);

  return returnVal;
}

//generalized lass lineshape as implemented in MINT3 by Tim Evans. if F=R=1 and phiF=PhiR=0 this is equal to lass_Mint3.
EXEC_TARGET devcomplex<fptype> glass_MINT3 (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital             = indices[4];
  fptype frFactor               = 1;
  fptype rMass2                 = Mpair*Mpair;
  
  fptype a = 2.07;
  fptype r = 3.32;
  fptype F = 1.0;
  fptype R = 1.0;
  fptype phiF = 0.0;
  fptype phiR = 0.0;

  fptype mpsq = (m1+m2)*(m1+m2);
  fptype mmsq = (m1-m2)*(m1-m2);
  fptype num  = (rMass2 - mpsq)*(rMass2 - mmsq);
  fptype num2  = (resmass*resmass - mpsq)*(resmass*resmass - mmsq);
  fptype pABSq = num/(4*rMass2);
  fptype prSq = FABS(num2/(4*resmass*resmass));

  fptype pratio = SQRT(pABSq/prSq);

  fptype pratio_to_2Jplus1 = 1;

  for(int i=0; i < 2*orbital+1; i++){
    pratio_to_2Jplus1 *= pratio;
  }

  fptype mratio = resmass/Mpair;
  fptype r2 = meson_radius*meson_radius;
  fptype thisFR = BL_PRIME(pABSq*r2, prSq*r2, orbital);
  fptype GofM = reswidth * pratio_to_2Jplus1 *mratio * thisFR;


  fptype y = 2.0 * a*SQRT(pABSq);
  fptype x = 2.0 + a * r * pABSq;
  fptype scattphase = phiF + atan(y/x);
  fptype resphase = phiR + atan(resmass*GofM/(resmass*resmass - rMass2)) ;
  fptype rho = 1.0 / SQRT(pABSq/rMass2);
  devcomplex<fptype> returnVal = (F * SIN(scattphase) * devcomplex<fptype>(COS(scattphase),SIN(scattphase)) 
                               + R * SIN(resphase) * devcomplex<fptype>(COS(resphase + 2 * scattphase),SIN(resphase + 2 * scattphase)))
                               * rho;

  // printf("GLass: %.5g %.5g %.5g\n",returnVal.real, returnVal.imag, GofM);
  return returnVal;
}




EXEC_TARGET devcomplex<fptype> aSqrtTerm(const fptype& m0, const fptype& m){
  fptype a2 = 1 - (2*m0/m)*(2*m0/m);
  devcomplex<fptype> returnVal = a2>0 ? devcomplex<fptype>(SQRT(a2),0) : devcomplex<fptype>(0, SQRT(-a2));
  return returnVal;
}

EXEC_TARGET devcomplex<fptype> Flatte_MINT (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital             = indices[4];
  fptype frFactor               = 1;
  fptype rMass2                 = Mpair*Mpair;
  
  // As far as I understand, this is only valid for the f980
  fptype gPi          = .165;
  fptype gK_by_gPi    = 4.21;
  fptype gK           = gPi * gK_by_gPi;
  fptype mPi0         = .1349766;
  fptype mPiPlus      = .13957018;
  fptype mKPlus       = .493677;
  fptype mK0          = .497648;

  fptype mpsq = (m1+m2)*(m1+m2);
  fptype mmsq = (m1-m2)*(m1-m2);
  fptype num  = (rMass2 - mpsq)*(rMass2 - mmsq);
  fptype num2  = (resmass*resmass - mpsq)*(resmass*resmass - mmsq);
  fptype pABSq = num/(4*rMass2);
  fptype prSq = FABS(num2/(4*resmass*resmass));

  devcomplex<fptype> Gpipi = (1./3.) * aSqrtTerm(mPi0, Mpair) + (2./3.) * aSqrtTerm(mPiPlus, Mpair);
  devcomplex<fptype> GKK = (1./2.) * aSqrtTerm(mK0, Mpair) + (1./2.) * aSqrtTerm(mKPlus, Mpair);
  devcomplex<fptype> FlatteWidth = gPi * Gpipi + gK * GKK;
  // printf("%.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Gpipi.real, Gpipi.imag, GKK.real, GKK.imag, FlatteWidth.real, FlatteWidth.imag, Mpair, pABSq);

  frFactor = BL(pABSq*meson_radius*meson_radius, orbital);
  devcomplex<fptype> BW = SQRT(frFactor) / devcomplex<fptype>(resmass*resmass - rMass2,0) - devcomplex<fptype>(0,1) * resmass * FlatteWidth;
  return BW;
}


EXEC_TARGET devcomplex<fptype> nonres_DP (fptype m12, fptype m1, fptype m2, unsigned int* indices) {
  return devcomplex<fptype>(1, 0); 
}

MEM_DEVICE resonance_function_ptr ptr_to_BW_DP4 = BW;
MEM_DEVICE resonance_function_ptr ptr_to_lass = lass_MINT;
MEM_DEVICE resonance_function_ptr ptr_to_lass3 = lass_MINT3;
MEM_DEVICE resonance_function_ptr ptr_to_glass3 = glass_MINT3;
MEM_DEVICE resonance_function_ptr ptr_to_bugg_MINT = bugg_MINT;
MEM_DEVICE resonance_function_ptr ptr_to_SBW = SBW;
MEM_DEVICE resonance_function_ptr ptr_to_NONRES_DP = nonres_DP;
MEM_DEVICE resonance_function_ptr ptr_to_Flatte = Flatte_MINT;

Lineshape::Lineshape (string name,
						Variable* mass, 
						Variable* width, 
						unsigned int L, 
						unsigned int Mpair,
            LS kind, 
            FF FormFac)
  : GooPdf(0, name), _mass(mass), _width(width), _L(L), _Mpair(Mpair), _kind(kind), _FormFac(FormFac)
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
  pindices.push_back(Mpair); 
  pindices.push_back(enum_to_underlying(FormFac)); 

  switch(kind){
    
    case LS::BW:
      GET_FUNCTION_ADDR(ptr_to_BW_DP4);
      break;
    
    case LS::Lass:
      GET_FUNCTION_ADDR(ptr_to_lass);
      break;

    case LS::Bugg:
      GET_FUNCTION_ADDR(ptr_to_bugg_MINT);
      break;

    case LS::SBW:
      GET_FUNCTION_ADDR(ptr_to_SBW);
      break;

    case LS::Flatte:
      GET_FUNCTION_ADDR(ptr_to_Flatte);
      break;
    
    default:
    fprintf(stderr,"It seems that the requested lineshape is not implemented yet. Check LineshapesPdf.cu");
    exit(0);
  }

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

Amplitude::Amplitude(std::string uniqueDecayStr, Variable* ar, Variable* ai, std::vector<Lineshape*> LS, std::vector<SpinFactor*> SF, unsigned int nPerm)
  : _uniqueDecayStr(uniqueDecayStr),
    _ar(ar),
    _ai(ai),
    _LS(LS),
    _SF(SF),
    _nPerm(nPerm)
{}

bool Amplitude::operator==(const Amplitude& A) const{
  return _uniqueDecayStr==A._uniqueDecayStr 
                and _ar == A._ar 
                and _ai == A._ai 
                and _LS == A._LS 
                and _SF == A._SF 
                and _nPerm == A._nPerm;
}



