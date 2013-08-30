#ifndef GAUSSIAN_THRUST_FUNCTOR_HH
#define GAUSSIAN_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class GaussianThrustFunctor : public ThrustPdfFunctor {
public:
  GaussianThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
