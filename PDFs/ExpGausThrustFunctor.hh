#ifndef EXPGAUS_THRUST_FUNCTOR_HH
#define EXPGAUS_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class ExpGausThrustFunctor : public EngineCore {
public:
  ExpGausThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 

private:

};

#endif
