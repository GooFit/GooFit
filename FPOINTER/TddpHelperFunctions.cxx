__device__ fptype twoBodyCMmom (double rMassSq, fptype d1m, fptype d2m) {
  // For A -> B + C, calculate momentum of B and C in rest frame of A. 
  // PDG 38.16.

  fptype kin1 = 1 - POW(d1m+d2m, 2) / rMassSq;
  if (kin1 >= 0) kin1 = SQRT(kin1);
  else kin1 = 1;
  fptype kin2 = 1 - POW(d1m-d2m, 2) / rMassSq;
  if (kin2 >= 0) kin2 = SQRT(kin2);
  else kin2 = 1; 

  return 0.5*SQRT(rMassSq)*kin1*kin2; 
}


__device__ fptype dampingFactorSquare (fptype cmmom, int spin, fptype mRadius) {
  fptype square = mRadius*mRadius*cmmom*cmmom;
  fptype dfsq = 1 + square; // This accounts for spin 1
  if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.   

  // Spin 3 and up not accounted for. 
  return dfsq; 
}

__device__ bool inDalitz (fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3) {
  if (m12 < POW(dm1 + dm2, 2)) return false; // This m12 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m12 > POW(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 
  
  // Calculate energies of 1 and 3 particles in m12 rest frame. 
  fptype e1star = 0.5 * (m12 - dm2*dm2 + dm1*dm1) / SQRT(m12); 
  fptype e3star = 0.5 * (bigM*bigM - m12 - dm3*dm3) / SQRT(m12); 

  // Bounds for m13 at this value of m12.
  fptype minimum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) + SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 < minimum) return false;
  fptype maximum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) - SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 > maximum) return false;

  return true; 
}

__device__ fptype spinFactor (unsigned int spin, fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, fptype m12, fptype m13, fptype m23, unsigned int cyclic_index) {
  if (0 == spin) return 1; // Should not cause branching since every thread evaluates the same resonance at the same time. 
  /*
  // Copied from BdkDMixDalitzAmp
   
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass)); 
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug3Mass)); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 
    
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m12 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m23 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 

  // The above, collapsed into single tests where possible. 
  fptype _mA = (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass);
  fptype _mB = (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 

  fptype _mAC = (PAIR_23 == cyclic_index ? m13 : m23);
  fptype _mBC = (PAIR_12 == cyclic_index ? m13 : m12);
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 
  */

  // Copied from EvtDalitzReso, with assumption that pairAng convention matches pipipi0 from EvtD0mixDalitz.
  // Again, all threads should get the same branch. 
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m12 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 

  fptype massFactor = 1.0/_mAB;
  fptype sFactor = -1; 
  sFactor *= ((_mBC - _mAC) + (massFactor*(motherMass*motherMass - _mC*_mC)*(_mA*_mA-_mB*_mB)));
  if (2 == spin) {
    sFactor *= sFactor; 
    fptype extraterm = ((_mAB-(2*motherMass*motherMass)-(2*_mC*_mC))+massFactor*pow((motherMass*motherMass-_mC*_mC),2));
    extraterm *= ((_mAB-(2*_mA*_mA)-(2*_mB*_mB))+massFactor*pow((_mA*_mA-_mB*_mB),2));
    extraterm /= 3;
    sFactor -= extraterm;
  }
  return sFactor; 
}

__device__ devcomplex<fptype> plainBW (fptype m12, fptype m13, fptype m23, 
				       fptype resmass, fptype reswidth, fptype meson_radius, 
				       fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, 
				       unsigned int spin, unsigned int cyclic_index) {
  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor = 1;

  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius); 
  }
  
  // RBW evaluation
  fptype A = (resmass - rMassSq); 
  fptype B = resmass*reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor / SQRT(rMassSq);
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> ret(A*C, B*C); // Dropping F_D=1

  ret *= SQRT(frFactor); 
  fptype spinF = spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index); 

  //if ((internalDebug1 == threadIdx.x) && (internalDebug2 == blockIdx.x))
    //printf("Breit-Wigner: %f %f %f | %f %f %f\n", A, B, C, resmass, rMassSq, reswidth);
    //printf("Breit-Wigner: %i %i %f %f | (%f %f) %f\n", spin, cyclic_index, m12, m13, ret.real, ret.imag, spinF);

  ret *= spinF; 
  return ret; 
}

__device__ devcomplex<fptype> gaussian (fptype m12, fptype m13, fptype m23, fptype resmass, fptype reswidth, unsigned int cyclic_index) {
  // Notice sqrt - this function uses mass, not mass-squared like the other resonance types. 
  fptype massToUse = SQRT(PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  //reswidth *= root2; 
  //reswidth = 1.0 / reswidth; 
  massToUse -= resmass;
  massToUse /= reswidth;
  massToUse *= massToUse;
  fptype ret = EXP(-0.5*massToUse); 

  // Ignore factor 1/sqrt(2pi). 
  ret /= reswidth;
  //ret *= invRootPi;

  //if ((internalDebug1 == threadIdx.x) && (internalDebug2 == blockIdx.x))
  //printf("Gaussian: %f %f %f | %f %f %f\n", m12, m13, m23, resmass, reswidth, ret); 

  return devcomplex<fptype>(ret, 0); 
}

__device__ fptype hFun (double s, double daug2Mass, double daug3Mass) {
  //Last helper function
  const fptype _pi = 3.14159265359;
  double sm   = daug2Mass + daug3Mass;
  double SQRTs = sqrt(s);
  double k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);

  double val = ((2/_pi) * (k_s/SQRTs) * log( (SQRTs + 2*k_s)/(sm)));

  return val;
}


__device__ fptype dh_dsFun (double s, double daug2Mass, double daug3Mass) {
  //Yet another helper function
  const fptype _pi = 3.14159265359;
  double k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);
  
  double val = (hFun(s, daug2Mass, daug3Mass) * (1.0/(8.0*pow(k_s, 2)) - 1.0/(2.0 * s)) + 1.0/(2.0* _pi*s));
  return val;
}


__device__ fptype dFun (double s, double daug2Mass, double daug3Mass) {
  // Helper function used in Gronau-Sakurai
  // static double _pi = TMath::Pi();
  const fptype _pi = 3.14159265359;
  double sm   = daug2Mass + daug3Mass;
  double sm24 = sm*sm/4.0;
  double m    = sqrt(s);
  double k_m2 = twoBodyCMmom(s, daug2Mass, daug3Mass);
 
  double val = 3.0/_pi * sm24/pow(k_m2, 2) * log((m + 2*k_m2)/sm) + m/(2*_pi*k_m2) - sm24*m/(_pi * pow(k_m2, 3));
 
  return val;
}

__device__ fptype fsFun (double s, double m2, double gam, double daug2Mass, double daug3Mass) {
  // Another G-S helper function
   
  double k_s   = twoBodyCMmom(s,  daug2Mass, daug3Mass);
  double k_Am2 = twoBodyCMmom(m2, daug2Mass, daug3Mass);
   
  double f     = gam * m2 / POW(k_Am2, 3);
  f           *= (POW(k_s, 2) * (hFun(s, daug2Mass, daug3Mass) - hFun(m2, daug2Mass, daug3Mass)) + (m2 - s) * pow(k_Am2, 2) * dh_dsFun(m2, daug2Mass, daug3Mass));
 
  return f;
}

__device__ devcomplex<fptype> gouSak (fptype m12, fptype m13, fptype m23, 
				       fptype resmass, fptype reswidth, fptype meson_radius, 
				       fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, 
				       unsigned int spin, unsigned int cyclic_index) {
  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor = 1;

  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius); 
  }
  
  // Implement Gro-Sak:

  fptype D = (1.0 + dFun(resmass, daug2Mass, daug3Mass) * reswidth/SQRT(resmass));
  fptype E = resmass - rMassSq + fsFun(rMassSq, resmass, reswidth, daug2Mass, daug3Mass);
  fptype F = SQRT(resmass) * reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor;

  D       /= (E*E + F*F);
  devcomplex<fptype> retur(D*E, D*F); // Dropping F_D=1
  retur *= SQRT(frFactor);
  retur *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);

  return retur; 
}


__device__ devcomplex<fptype> getResonanceAmplitude (fptype m12, fptype m13, fptype resmass, fptype reswidth, unsigned int spin, unsigned int cyclic_index, unsigned int eval_type, unsigned int* indices) {
  fptype motherMass   = functorConstants[indices[1] + 0]; 
  fptype daug1Mass    = functorConstants[indices[1] + 1]; 
  fptype daug2Mass    = functorConstants[indices[1] + 2]; 
  fptype daug3Mass    = functorConstants[indices[1] + 3]; 
  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return devcomplex<fptype>(0, 0); 
  fptype meson_radius = functorConstants[indices[1] + 4]; 

  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 

  // This switch evaluates identically for each thread, since we're going through the resonances
  // in sequence. So there should not be any branching because of it. 
  switch (eval_type) {
  case RBW:    return plainBW(m12, m13, m23, resmass, reswidth, meson_radius, motherMass, daug1Mass, daug2Mass, daug3Mass, spin, cyclic_index); 
  case GOU_SAK: return  gouSak(m12, m13, m23, resmass, reswidth, meson_radius, motherMass, daug1Mass, daug2Mass, daug3Mass, spin, cyclic_index);
  case GAUSSIAN: return gaussian(m12, m13, m23, resmass, reswidth, cyclic_index); 
  default: 
  case NONRES: return devcomplex<fptype>(1, 0); 
  }  
}

__device__ void getAmplitudeCoefficients (devcomplex<fptype> a1, devcomplex<fptype> a2, fptype& a1sq, fptype& a2sq, fptype& a1a2real, fptype& a1a2imag) {
  // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
  a1sq = a1.abs2();
  a2sq = a2.abs2();
  a1 *= conj(a2);
  a1a2real = a1.real;
  a1a2imag = a1.imag; 
}
