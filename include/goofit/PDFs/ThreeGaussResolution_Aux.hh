#ifndef THREEGAUSS_RESOLUTION_HH
#define THREEGAUSS_RESOLUTION_HH

#include "MixingTimeResolution_Aux.hh"

class ThreeGaussResolution : public MixingTimeResolution {
public: 
  ThreeGaussResolution (Variable* cf, Variable* tf, Variable* cb, Variable* cs, Variable* tb, Variable* ts, Variable* ob, Variable* os); 
  ~ThreeGaussResolution ();

  virtual fptype normalisation (fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const;
  virtual void createParameters (std::vector<unsigned int>& pindices, PdfBase* dis); 

private:
  Variable* coreFraction;
  Variable* tailFraction;
  Variable* coreBias;
  Variable* coreScaleFactor;
  Variable* tailBias;
  Variable* tailScaleFactor;
  Variable* outBias;
  Variable* outScaleFactor;
}; 

#endif 
