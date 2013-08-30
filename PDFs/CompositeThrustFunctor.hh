#ifndef COMPOSITE_THRUST_FUNCTOR_HH
#define COMPOSITE_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

// Composites of arbitrary functions, ie f(x) = h(g(x)) 
// for any h and g. In principle we should allow multi-
// dimensional compositing, eg f(x, y) = i(g(x, y), h(x, y)).
// Not implemented yet. 

class CompositeThrustFunctor : public ThrustPdfFunctor {
public:
  CompositeThrustFunctor (std::string n, FunctorBase* core, FunctorBase* shell); // Where 'core' corresponds to 'g' and 'shell' to 'h'. 
  __host__ virtual fptype normalise () const;

private:

};

#endif
