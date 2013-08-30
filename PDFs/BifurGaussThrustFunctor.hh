#ifndef BIFURGAUSSIAN_THRUST_FUNCTOR_HH
#define BIFURGAUSSIAN_THRUST_FUNCTOR_HH

#include "EngineCore.hh"

class BifurGaussThrustFunctor : public EngineCore {
  public:
    BifurGaussThrustFunctor (std::string n, Variable *_x, Variable* m, Variable* sL, Variable* sR);
    __host__ fptype integrate(fptype lo, fptype hi) const;
    //__host__ virtual bool hasAnalyticIntegral () const {return true;}

  private:
 
};

#endif
