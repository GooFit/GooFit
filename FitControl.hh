#ifndef FITCONTROL_HH
#define FITCONTROL_HH

#include <string> 
#include "GlobalCudaDefines.hh" 
#include "Variable.hh" 
#include <vector> 

class FunctorBase; 

class FitControl {
public:
  FitControl (bool bin, std::string mn); 
  ~FitControl (); 

  inline bool binnedFit () const {return binned;}
  inline bool binErrors () const {return errorsOnBins;} 
  inline bool metricIsPdf () const {return !errorsOnBins;} 
  inline std::string getMetric () const {return metricName;} 
  inline FunctorBase* getOwner () const {return owner;} 
  void setOwner (FunctorBase* dat);

protected: 
  bool errorsOnBins; 

private:
  bool binned; 
  std::string metricName; 
  FunctorBase* owner; 
}; 

class UnbinnedNllFit : public FitControl {
public:
  UnbinnedNllFit ();
};

class BinnedNllFit : public FitControl {
public:
  BinnedNllFit ();
};

class BinnedErrorFit : public FitControl {
public:
  BinnedErrorFit (); 
};

class BinnedChisqFit : public FitControl {
public:
  BinnedChisqFit (); 
};



#endif
