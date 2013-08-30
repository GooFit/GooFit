#ifndef SCALEDGAUSSIAN_THRUST_FUNCTOR_HH
#define SCALEDGAUSSIAN_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class ScaledGaussianThrustFunctor : public ThrustPdfFunctor {
public:
  ScaledGaussianThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* d, Variable* e); 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 



private:

};

#endif
