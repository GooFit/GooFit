#ifndef LINESHAPES_PDF_HH
#define LINESHAPES_PDF_HH

#include "GooPdf.hh" 
#include "devcomplex.hh" 
#include "ResonancePdf.hh"

//typedef devcomplex<fptype> (*resonance_function_ptr) (fptype, fptype, fptype, unsigned int*); 
typedef fptype (*spin_function_ptr) (fptype*, unsigned int*); 
EXEC_TARGET void get4Vecs (fptype* Vecs, const unsigned int& constants, const fptype& m12, const fptype& m34, const fptype& cos12, const fptype& cos34, const fptype& phi);
EXEC_TARGET fptype getmass(const unsigned int& pair, fptype& d1, fptype& d2, const fptype* vecs, const fptype& m1, const fptype& m2, const fptype& m3, const fptype& m4);

class __align__(16) gpuLVec
{
  private:
    fptype X;
    fptype Y;
    fptype Z;
    fptype E;
  public:
    EXEC_TARGET gpuLVec(fptype x, fptype y, fptype z, fptype e);
    EXEC_TARGET fptype getX() const {return X;}
    EXEC_TARGET fptype getY() const {return Y;}
    EXEC_TARGET fptype getZ() const {return Z;}
    EXEC_TARGET fptype getE() const {return E;}
    EXEC_TARGET fptype Dot(const gpuLVec& rhs) const;
    EXEC_TARGET gpuLVec& operator+=(const gpuLVec& rhs);
    EXEC_TARGET gpuLVec& operator-=(const gpuLVec& rhs);
};

EXEC_TARGET gpuLVec operator+(gpuLVec lhs, const gpuLVec& rhs);
EXEC_TARGET gpuLVec operator-(gpuLVec lhs, const gpuLVec& rhs);


EXEC_TARGET fptype LeviCevita(const gpuLVec& p1, const gpuLVec& p2, const gpuLVec& p3, const gpuLVec& p4);
class SpinFactor;

class Lineshape : public GooPdf {
  // Service class intended to hold parametrisations of
  // resonances on Dalitz plots. Don't try to use this
  // as a standalone PDF! It should only be used as a
  // component in one of the friend classes. It extends
  // GooPdf so as to take advantage of the 
  // infrastructure, but will crash if used on its own. 

  friend class DPPdf; 
public:
  // Constructor for regular BW 
  Lineshape (string name,
			  Variable* mass, 
			  Variable* width, 
			  unsigned int L, 
			  unsigned int cyc); 

  // Nonresonant constructor
  Lineshape (string name);  
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}

private:
  unsigned int _mother_pdg;
  /*
  Variable* mass;
  Variable* width;
  unsigned int spin;
  unsigned int cyclic_index;
  unsigned int eval_type;
  unsigned int resonance_type; 
  */ 
};

class Amplitude {
  friend class DPPdf;

public:
  Amplitude(std::string uniqueDecayStr, Variable* ar, Variable* ai, std::map<std::string, Lineshape*> LS, std::map<std::string, SpinFactor*> SF);
  
private:
  std::string _uniqueDecayStr;
  Variable* _ar;
  Variable* _ai;
  std::map<std::string, SpinFactor*> _SF;
  std::map<std::string, Lineshape*> _LS;
};


class SpinFactor : public GooPdf {
  friend class DPPdf;

public:
  SpinFactor(std::string name, unsigned int kind, unsigned int P0, unsigned int P1, unsigned int P2, unsigned int P3);
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}
  
private:
  unsigned int kind;
};

// class SpinCalculator : public thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype> >, public GooPdf {
// public:
//   // Used to create the cached SF values. 
//   SpinCalculator (std::string name, unsigned int kind); 
//   EXEC_TARGET devcomplex<fptype> operator () (thrust::tuple<int, fptype*, int> t) const ;
//   // void SetParameterIdx(const unsigned int &idx) {parameters  = idx;}
//   // void resolveMassIdx(std::map<unsigned int, unsigned int> massmap );
// private:
//   unsigned int parameters;
//   unsigned int mother;
//   unsigned int self;
//   unsigned int d1;
//   unsigned int d2;

// }; 



#endif
