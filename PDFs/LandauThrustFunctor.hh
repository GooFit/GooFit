#ifndef LANDAU_THRUST_FUNCTOR_HH
#define LANDAU_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class LandauThrustFunctor : public EngineCore {
public:
  LandauThrustFunctor (std::string n, Variable* _x, Variable* mpv, Variable* sigma); 

private:

};

#endif
