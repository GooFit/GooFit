#ifndef ARGUS_THRUST_FUNCTOR_HH
#define ARGUS_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class ArgusThrustFunctor : public ThrustPdfFunctor {
public:
  ArgusThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, bool upper, Variable* power = 0); 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 
  __host__ fptype integrate (fptype lo, fptype hi) const; 

private:

};

#endif
