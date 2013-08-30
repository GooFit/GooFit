#ifndef STEP_THRUST_FUNCTOR_HH
#define STEP_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class StepThrustFunctor : public EngineCore {
public:
  StepThrustFunctor (std::string n, Variable* _x, Variable* x0); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
