/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#ifndef LINESHAPES_PDF_HH
#define LINESHAPES_PDF_HH

#include "GooPdf.hh" 
#include "devcomplex.hh" 
#include "ResonancePdf.hh"
#include "SpinFactors.hh"

enum class LS {BW, Lass, Bugg, Flatte, SBW};
//PDG notation for FF
enum class FF : unsigned int {One = 0, BL, BL_Prime};

class Lineshape : public GooPdf {
  friend class DPPdf; 
  // Service class intended to hold parametrisations of
  // resonances on Dalitz plots. Don't try to use this
  // as a standalone PDF! It should only be used as a
  // component in one of the friend classes. It extends
  // GooPdf so as to take advantage of the 
  // infrastructure, but will crash if used on its own. 
  Variable* _mass;
  Variable* _width; 
  unsigned int _L; 
  unsigned int _Mpair;
  LS _kind;
  FF _FormFac;
public:
  Lineshape (string name,
        Variable* mass, 
        Variable* width, 
        unsigned int L, 
        unsigned int Mpair,
        LS kind = LS::BW,
        FF FormFac = FF::BL_Prime); 

  bool operator==(const Lineshape& L) const {return ( L.getName() == getName() and L._mass->value == _mass->value and L._width->value == _width->value
                                                      and L._L == _L and L._Mpair == _Mpair and L._kind == _kind and L._FormFac == _FormFac); }
  Lineshape (string name);
   
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}
};

class Amplitude {
  friend class DPPdf;

public:
  Amplitude(std::string uniqueDecayStr, Variable* ar, Variable* ai, std::vector<Lineshape*> LS, std::vector<SpinFactor*> SF, unsigned int nPerm = 1);
  // bool operator==(const Amplitude& A){return ( A._uniqueDecayStr == _uniqueDecayStr ); }
private:
  std::string _uniqueDecayStr;
  Variable* _ar;
  Variable* _ai;
  std::vector<SpinFactor*> _SF;
  std::vector<Lineshape*> _LS;
  unsigned int _nPerm;
};

#endif
