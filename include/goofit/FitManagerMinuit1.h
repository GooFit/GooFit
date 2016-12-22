#ifndef FITMANAGER_MINUIT1_HH
#define FITMANAGER_MINUIT1_HH

#include "goofit/rootstuff/TMinuit.h" 
extern PdfBase* pdfPointer; 
extern int numPars; 

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag); 

class FitManager { 
public:
  FitManager (PdfBase* dat);
  ~FitManager ();
  void setMaxCalls (double mxc) {overrideCallLimit = mxc;}
  void useHesse (bool use=true) {_useHesse = use;}
  void useMinos (bool use=true) {_useMinos = use;}
  void useImprove (bool use=true) {_useImprove = use;}
  void setupMinuit ();
  void runMigrad (); 
  void fit (); 
  TMinuit* getMinuitObject () {return minuit;} 
  void getMinuitValues () const;
  void getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) const;
  TMinuit* minuit; 
private:
  double overrideCallLimit; 
  bool _useHesse;
  bool _useMinos;
  bool _useImprove;
};

#endif 
