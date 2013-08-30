#ifndef KINLIMITBW_THRUST_FUNCTOR_HH
#define KINLIMITBW_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class KinLimitBWThrustFunctor : public EngineCore {

public:
  KinLimitBWThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s); 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 
  __host__ void setMasses (fptype bigM, fptype smallM); 

private:

};

#endif
