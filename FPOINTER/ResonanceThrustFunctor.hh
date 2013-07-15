#ifndef RESONANCE_THRUST_FUNCTOR_HH
#define RESONANCE_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 
#include "devcomplex.hh" 
typedef devcomplex<fptype> (*resonance_function_ptr) (fptype, fptype, fptype, unsigned int*); 

class ResonanceThrustFunctor : public ThrustPdfFunctor {
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
