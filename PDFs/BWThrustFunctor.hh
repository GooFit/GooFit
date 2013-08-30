#ifndef BW_THRUST_FUNCTOR_HH
#define BW_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class BWThrustFunctor : public ThrustPdfFunctor {

public:
  BWThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s); 
private:

};

#endif
