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
			  unsigned int cyc,
        bool useMINTBW); 

  // Nonresonant constructor
  Lineshape (string name);  
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}

private:
  /*
  unsigned int _mother_pdg;
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
  Amplitude(std::string uniqueDecayStr, Variable* ar, Variable* ai, std::map<std::string, Lineshape*> LS, std::map<std::string, SpinFactor*> SF, unsigned int nPerm = 1);
  
private:
  std::string _uniqueDecayStr;
  Variable* _ar;
  Variable* _ai;
  std::map<std::string, SpinFactor*> _SF;
  std::map<std::string, Lineshape*> _LS;
  unsigned int _nPerm;
};

#endif
