#ifndef DALITZVETO_THRUST_FUNCTOR_HH
#define DALITZVETO_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 
#include "TddpThrustFunctor.hh"

struct VetoInfo {
  DaughterPair cyclic_index; 
  Variable* minimum;
  Variable* maximum; 
};

class DalitzVetoThrustFunctor : public EngineCore {
public:
  __host__ DalitzVetoThrustFunctor (std::string n,  Variable* _x, Variable* _y, Variable* motherM, Variable* d1m, Variable* d2m, Variable* d3m, vector<VetoInfo*> vetos);

private:

};

#endif
