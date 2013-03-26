#ifndef NOVOSIBIRSK_THRUST_FUNCTOR_HH
#define NOVOSIBIRSK_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class NovosibirskThrustFunctor : public ThrustPdfFunctor {
public:
  NovosibirskThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 

private:

};

#endif
