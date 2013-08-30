#ifndef LANDAU_THRUST_FUNCTOR_HH
#define LANDAU_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class LandauThrustFunctor : public ThrustPdfFunctor {
public:
  LandauThrustFunctor (std::string n, Variable* _x, Variable* mpv, Variable* sigma); 

private:

};

#endif
