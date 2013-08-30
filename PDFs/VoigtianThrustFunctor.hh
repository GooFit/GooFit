#ifndef VOIGTIAN_THRUST_FUNCTOR_HH
#define VOIGTIAN_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class VoigtianThrustFunctor : public ThrustPdfFunctor {
public:
  VoigtianThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* w); 

private:

};

#endif
