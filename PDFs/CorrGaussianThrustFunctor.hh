#ifndef CORRGAUSSIAN_THRUST_FUNCTOR_HH
#define CORRGAUSSIAN_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 

class CorrGaussianThrustFunctor : public ThrustPdfFunctor {
public:
  CorrGaussianThrustFunctor (std::string n, Variable* _x, Variable* _y, Variable* mean1, Variable* sigma1, Variable* mean2, Variable* sigma2, Variable* correlation); 


private:

};

#endif
