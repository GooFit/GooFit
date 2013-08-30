#ifndef BINTRANSFORM_THRUST_FUNCTOR_HH
#define BINTRANSFORM_THRUST_FUNCTOR_HH

#include "EngineCore.hh" 

// Transforms ND coordinates into a single bin number. 
class BinTransformThrustFunctor : public EngineCore {
public:
  BinTransformThrustFunctor (std::string n, vector<Variable*> obses, vector<fptype> limits, vector<fptype> binSizes, vector<int> numBins); 

private:

};

#endif
