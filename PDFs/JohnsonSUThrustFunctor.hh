#ifndef JOHNSONSU_THRUST_FUNCTOR_HH
#define JOHNSONSU_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class JohnsonSUThrustFunctor : public ThrustPdfFunctor {
public:
  JohnsonSUThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* g, Variable* d); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
