#ifndef ADD_THRUST_FUNCTOR_HH
#define ADD_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class AddThrustFunctor : public ThrustPdfFunctor {
public:

  AddThrustFunctor (std::string n, std::vector<Variable*> weights, std::vector<FunctorBase*> comps); 
  AddThrustFunctor (std::string n, Variable* frac1, FunctorBase* func1, FunctorBase* func2); 
  __host__ virtual fptype normalise () const;
  __host__ virtual bool hasAnalyticIntegral () const {return false;}

protected:
  __host__ virtual double sumOfNll (int numVars) const;

private:
  bool extended; 
};

#endif
