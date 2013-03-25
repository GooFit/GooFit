#ifndef CRYSTALBALL_THRUST_FUNCTOR_HH
#define CRYSTALBALL_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class CrystalBallThrustFunctor : public ThrustPdfFunctor {
public:
  CrystalBallThrustFunctor (std::string n, Variable* _x, Variable* m, Variable* s, Variable* a, Variable* power = 0); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  //__host__ virtual bool hasAnalyticIntegral () const {return true;} 


private:

};

#endif
