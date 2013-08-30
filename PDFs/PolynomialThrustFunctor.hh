#ifndef POLYNOMIAL_THRUST_FUNCTOR_HH
#define POLYNOMIAL_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class PolynomialThrustFunctor : public EngineCore {
public:
  PolynomialThrustFunctor (std::string n, Variable* _x, std::vector<Variable*> weights, Variable* x0 = 0, unsigned int lowestDegree = 0); 
  PolynomialThrustFunctor (string n, vector<Variable*> obses, vector<Variable*> coeffs, vector<Variable*> offsets, unsigned int maxDegree); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  //__host__ virtual bool hasAnalyticIntegral () const {return (1 == observables.size());} 
  __host__ fptype getCoefficient (int coef) const;


private:
  Variable* center; 
};

#endif
