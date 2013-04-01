#ifndef EXPGAUS_THRUST_FUNCTOR_HH
#define EXPGAUS_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class ExpGausThrustFunctor : public ThrustPdfFunctor {
public:
  ExpGausThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 

private:

};

#endif
