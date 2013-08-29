#ifndef RESONANCE_THRUST_FUNCTOR_HH
#define RESONANCE_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 
#include "devcomplex.hh" 
typedef devcomplex<fptype> (*resonance_function_ptr) (fptype, fptype, fptype, unsigned int*); 

class ResonanceThrustFunctor : public ThrustPdfFunctor {
  // Service class intended to hold parametrisations of
  // resonances on Dalitz plots. Don't try to use this
  // as a standalone PDF! It should only be used as a
  // component in one of the friend classes. It extends
  // ThrustPdfFunctor so as to take advantage of the 
  // infrastructure, but will crash if used on its own. 

  friend class TddpThrustFunctor;
  friend class DalitzPlotThrustFunctor; 
  friend class IncoherentSumThrustFunctor; 
public:
  // Constructor for regular BW 
  ResonanceThrustFunctor (string name, 
			  Variable* ar, 
			  Variable* ai, 
			  Variable* mass, 
			  Variable* width, 
			  unsigned int sp, 
			  unsigned int cyc); 

  // Gounaris-Sakurai
  ResonanceThrustFunctor (string name, 
			  Variable* ar, 
			  Variable* ai, 
			  unsigned int sp, 
			  Variable* mass, 
			  Variable* width, 
			  unsigned int cyc); 
 
  // LASS constructor
  ResonanceThrustFunctor (string name,
                          Variable* ar,
                          Variable* ai,
			  Variable* mass,
			  unsigned int sp,
                          Variable* width,
                          unsigned int cyc);
  

  // Nonresonant constructor
  ResonanceThrustFunctor (string name, 
			  Variable* ar, 
			  Variable* ai);  

  // Gaussian constructor
  ResonanceThrustFunctor (string name,
			  Variable* ar, 
			  Variable* ai,
			  Variable* mean, 
			  Variable* sigma,
			  unsigned int cyc);

private:
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}

  Variable* amp_real;
  Variable* amp_imag;
  /*
  Variable* mass;
  Variable* width;
  unsigned int spin;
  unsigned int cyclic_index;
  unsigned int eval_type;
  unsigned int resonance_type; 
  */ 
};

#endif
