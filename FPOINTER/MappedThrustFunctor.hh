#ifndef MAPPED_THRUST_FUNCTOR_HH
#define MAPPED_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class MappedThrustFunctor : public ThrustPdfFunctor {
public:
  MappedThrustFunctor (std::string n, ThrustPdfFunctor* m, vector<ThrustPdfFunctor*>& t); 
  // Map function m must be custom written to correspond to order of function list t. 
  __host__ fptype normalise () const;
private:

};

#endif
