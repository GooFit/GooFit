#ifndef TRIGTHRESHOLD_THRUST_FUNCTOR_HH
#define TRIGTHRESHOLD_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class TrigThresholdThrustFunctor : public ThrustPdfFunctor {
public:
  TrigThresholdThrustFunctor (std::string n, Variable* _x, Variable* thresh, Variable* trigConst, Variable* linConst, bool upper = true); 
  TrigThresholdThrustFunctor (std::string n, Variable* _x, Variable* _y, Variable* thresh, Variable* trigConst, Variable* linConst, Variable* massConstant, bool upper);

private:

};

#endif
