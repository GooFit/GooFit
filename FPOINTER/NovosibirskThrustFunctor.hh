#ifndef NOVOSIBIRSK_THRUST_FUNCTOR_HH
#define NOVOSIBIRSK_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class NovosibirskThrustFunctor : public ThrustPdfFunctor {
public:
  NovosibirskThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 
  /// __host__ fptype integrate (fptype lo, fptype hi) const; 
  /// __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
