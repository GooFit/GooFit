/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

TODO:
- reorganize this file into multiple files so lineshapes and spinfactors are seperated. Also all helper functions should go into seperate file.

This file includes some lineshapes and spinfactors.
Also right now it is the home to some helper functions needed and an implementation of a simple 4-vec class that works on the GPU
*/

#include "LineshapesPdf.hh" 

EXEC_TARGET fptype BL_PRIME (fptype z2, fptype z02, int L) {
  if (1 == L) return(1+z02)/(1+z2);
  else if (2 == L) return ( z02*z02 + 3*z02 + 9 ) / ( z2*z2 + 3*z2 + 9 )  ;
  else{
    printf("ERROR! Oribtal > 2 not supported!\n");
    return 0;
  }
  // Spin 3 and up not accounted for. 
}

EXEC_TARGET fptype BL(fptype z2, int L) {
  if( 1==L) return 2*z2/(1+z2);
  else if( 1==L) return (13*z2*z2 )/ ( z2*z2 + 3*z2 + 9 )  ;
  else{
    printf("ERROR! Oribtal > 2 not supported!\n");
    return 0;
  }
  // Spin 3 and up not accounted for. 
}

EXEC_TARGET devcomplex<fptype> BW_DP (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
  fptype meson_radius           = functorConstants[indices[1]+0];
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int orbital          = indices[4];
  unsigned int FF               = indices[6];

  fptype frFactor = 1;

  fptype rMass2 = Mpair * Mpair;
  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame 
  fptype q = twoBodyCMmom(rMass2, m1, m2);
  fptype q0 = twoBodyCMmom(resmass, m1, m2);

  fptype r2 = meson_radius*meson_radius;
  fptype q2 = q*q;

  if (0 != orbital and 0 != FF) {
    frFactor =  (FF==1? BL(q2*r2, orbital) : BL_PRIME(q2*r2, q0*q0*r2, orbital));
  }  
  
  fptype FF_Gamma = (2==FF ? frFactor : BL_PRIME(q2*r2, q0*q0*r2, orbital));
  // RBW evaluation
  fptype A = (resmass - rMass2); 
  fptype B = resmass*reswidth * POW(q / q0, 2.0*orbital + 1) * frFactor / SQRT(rMass2);
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> ret(A*C, B*C); // Dropping F_D=1

  ret *= SQRT(frFactor); 
  return ret; 
}

//This function is modeled after BW_BW::getVal() in BW_BW.cpp from the MINT package written by Jonas Rademacker. 
EXEC_TARGET devcomplex<fptype> BW_MINT (fptype Mpair, fptype m1, fptype m2, unsigned int* indices) {
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
  fptype prSqForGofM = FABS(num2/(4*mass*mass));
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
    frFactor =  (FF==1? BL(pABSq*r2, orbital) : BL_PRIME(pABSq*r2, prSqForGofM*r2, orbital));
  } 

  fptype GofM = width * pratio_to_2Jplus1 *mratio * thisFR * thisFR;

  fptype gamma = SQRT(mass*mass*(mass*mass + width*width));
  fptype k     = mass*width*gamma/SQRT(mass*mass+gamma);

  devcomplex<fptype> BW(mass*mass - mumsRecoMass2, mass*GofM);
  fptype den = (mass*mass - mumsRecoMass2) * (mass*mass - mumsRecoMass2) + mass * GofM * mass * GofM;

  devcomplex<fptype> ret = (SQRT(k) * frFactor)/den * BW;

  // printf("%.7g, %.7g, %.7g, %i, %.7g, %.7g\n", meson_radius, resmass, reswidth, orbital, pABSq, prSqForGofM);
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

//This function is an adaptation from the bugg lineshape implemented in the MINT package written by Jonas Rademacker. 
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
MEM_DEVICE resonance_function_ptr ptr_to_bugg_MINT = bugg_MINT;
MEM_DEVICE resonance_function_ptr ptr_to_NONRES_DP = nonres_DP;

Lineshape::Lineshape (string name,
						Variable* mass, 
						Variable* width, 
						unsigned int L, 
						unsigned int Pair,
            LS kind, 
            FF FormFac)
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
  pindices.push_back(Pair); 
  pindices.push_back(enum_to_underlying(FormFac)); 

  switch(kind){
    
    case LS::BW:
      GET_FUNCTION_ADDR(ptr_to_BW_DP);
      break;
    
    case LS::BW_MINT:
      GET_FUNCTION_ADDR(ptr_to_BW_MINT);
      break;

    case LS::Lass:
      GET_FUNCTION_ADDR(ptr_to_lass_DP);
      break;

    case LS::Bugg:
    GET_FUNCTION_ADDR(ptr_to_bugg_MINT);
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

Amplitude::Amplitude(std::string uniqueDecayStr, Variable* ar, Variable* ai, std::map<std::string, Lineshape*> LS, std::map<std::string, SpinFactor*> SF, unsigned int nPerm)
  : _uniqueDecayStr(uniqueDecayStr),
    _ar(ar),
    _ai(ai),
    _LS(LS),
    _SF(SF),
    _nPerm(nPerm)
{}



