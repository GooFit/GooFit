#ifndef PROD_THRUST_FUNCTOR_HH
#define PROD_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

class ProdThrustFunctor : public EngineCore {
public:

  ProdThrustFunctor (std::string n, std::vector<FunctorBase*> comps); 
  __host__ virtual fptype normalise () const;
  __host__ virtual bool hasAnalyticIntegral () const {return false;}

private:
  bool varOverlaps; // True if any components share an observable. 
};

#endif
