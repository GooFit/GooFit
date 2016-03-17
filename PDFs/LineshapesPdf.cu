#include "LineshapesPdf.hh" 


// EXEC_TARGET fptype twoBodyCMmom (double rMassSq, fptype d1m, fptype d2m)
// defined in ResonancePdf

//EXEC_TARGET fptype dampingFactorSquare (fptype cmmom, int spin, fptype mRadius) {
// defined in ResonancePdf

//EXEC_TARGET devcomplex<fptype> spin0 (fptype** Vecs, unsigned int* indices) {
EXEC_TARGET devcomplex<fptype> spin0 (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  return devcomplex<fptype>(1,0);
}

EXEC_TARGET devcomplex<fptype> BW_DP (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+1];
  fptype daug1Mass              = functorConstants[indices[1]+2];
  fptype daug2Mass              = functorConstants[indices[1]+3];
  fptype daug3Mass              = functorConstants[indices[1]+4];
  fptype meson_radius           = functorConstants[indices[1]+0];

  fptype resmass                = cudaArray[indices[3]];
  fptype reswidth               = cudaArray[indices[4]];
  unsigned int orbital          = indices[5];
  unsigned int cyclic_index     = indices[6]; 

//  fptype m12 = (P1[0]+P2[0]) * (P1[0]+P2[0]) + (P1[1]+P2[1]) * (P1[1]+P2[1]) + (P1[2]+P2[2]) * (P1[2]+P2[2]) + (P1[3]+P2[3]) * (P1[3]+P2[3]);
//  fptype m1 = P1[0] * P1[0] + P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3];
//  fptype m2 = P2[0] * P2[0] + P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3];
  fptype frFactor = 1;

  // printf("%f, %f, %i, %i\n", resmass, reswidth, orbital, cyclic_index);
  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(m12, daug1Mass, daug2Mass);
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, daug1Mass, daug2Mass);

  if (0 != orbital) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, orbital, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, orbital, meson_radius); 
  }  
  // RBW evaluation
  fptype A = (resmass - m12); 
  fptype B = resmass*reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*orbital + 1) * frFactor / SQRT(m12);
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> ret(A*C, B*C); // Dropping F_D=1

  fptype sf = spinFactor(orbital, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);
  ret *= SQRT(frFactor); 
  ret *=sf;
  printf("%f, %f, %f, %f\n",ret.real, ret.imag, m12, m13);
  return ret; 
}

EXEC_TARGET devcomplex<fptype> nonres_DP (fptype m12, fptype m1, fptype m2, unsigned int* indices) {
  return devcomplex<fptype>(1, 0); 
}

MEM_DEVICE resonance_function_ptr ptr_to_BW_DP = BW_DP;
MEM_DEVICE resonance_function_ptr ptr_to_NONRES_DP = nonres_DP;
MEM_DEVICE resonance_function_ptr ptr_to_spin0 = spin0;



Lineshape::Lineshape (string name,
            unsigned int mother_pdg, 
						Variable* mass, 
						Variable* width, 
						unsigned int sp, 
						unsigned int cyc) 
  : GooPdf(0, name),
    _mother_pdg(mother_pdg)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Making room for index of decay-related constants. Assumption:
  // These are mother mass and three daughter masses in that order.
  // They will be registered by the object that uses this resonance,
  // which will tell this object where to find them by calling setConstantIndex. 
  pindices.push_back(0); //index of mother

  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(sp);
  pindices.push_back(cyc); 

  GET_FUNCTION_ADDR(ptr_to_BW_DP);
  initialise(pindices); 
}



Lineshape::Lineshape (string name, unsigned int mother_pdg)
  : GooPdf(0, name),
    _mother_pdg(mother_pdg)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  GET_FUNCTION_ADDR(ptr_to_NONRES_DP);
  initialise(pindices); 
}

Amplitude::Amplitude(std::string uniqueDecayStr, Variable* ar, Variable* ai, std::map<std::string, Lineshape*> LS, std::map<std::string, SpinFactor*> SF )
  : _uniqueDecayStr(uniqueDecayStr),
    _ar(ar),
    _ai(ai),
    _LS(LS),
    _SF(SF)
{}



SpinFactor::SpinFactor (std::string name, unsigned int kind, unsigned int P0, unsigned int P1, unsigned int P2, unsigned int P3)
 : GooPdf(0,name){
  vector<unsigned int> pindices; 
  pindices.push_back(0); //dummy for index to constants.
  pindices.push_back(P0);
  pindices.push_back(P1);
  pindices.push_back(P2);
  pindices.push_back(P3);
  if (kind==0){
    GET_FUNCTION_ADDR(ptr_to_spin0);
  }else {
    std::cout << "No Spinfunction implemented for that kind." << std::endl;
    exit(0);
  }
  
  initialise(pindices);
}

// EXEC_TARGET devcomplex<fptype> SpinCalculator::operator () (thrust::tuple<int, fptype*, int> t) const{
//   return devcomplex<fptype>(1,0);
// }

// void SpinCalculator::resolveMassIdx(std::map<unsigned int, unsigned int> massmap ){
//   mother = massmap[mother];
//   self = massmap[self];
//   d1 = massmap[d1];
//   d2 = massmap[d2];
// }
