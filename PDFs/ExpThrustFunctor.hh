#ifndef EXP_THRUST_FUNCTOR_HH
#define EXP_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class ExpThrustFunctor : public EngineCore {
public:
  ExpThrustFunctor (std::string n, Variable* _x, Variable* alpha, Variable* offset = 0); 
  ExpThrustFunctor (std::string n, Variable* _x, std::vector<Variable*>& weights, Variable* offset = 0); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return (1 == host_indices[parameters]);} 



private:

};

#endif
