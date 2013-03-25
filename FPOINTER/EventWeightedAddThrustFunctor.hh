#ifndef EVENTWEIGHTEDADD_THRUST_FUNCTOR_HH
#define EVENTWEIGHTEDADD_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

// This class is just like AddThrustFunctor except that the
// event weights are properties of each event, not variables
// in the fit. 
class EventWeightedAddThrustFunctor : public ThrustPdfFunctor {
public:

  EventWeightedAddThrustFunctor (std::string n, std::vector<Variable*> weights, std::vector<FunctorBase*> comps); 
  __host__ virtual fptype normalise () const;
  __host__ virtual bool hasAnalyticIntegral () const {return false;}

protected:

private:

};

#endif
