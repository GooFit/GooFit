#ifndef PDFBUILDER_MINUIT1_HH
#define PDFBUILDER_MINUIT1_HH

#include "TMinuit.hh" 
extern FunctorBase* pdfPointer; 
extern int numPars; 
#ifdef OMP_ON
#pragma omp threadprivate (numPars)
#pragma omp threadprivate (pdfPointer)
#endif

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag); 

class PdfFunctor { 
public:
  PdfFunctor (FunctorBase* dat);
  ~PdfFunctor ();
  void setMaxCalls (double mxc) {overrideCallLimit = mxc;}
  void setupMinuit ();
  void runMigrad (); 
  void fit (); 
  TMinuit* getMinuitObject () {return minuit;} 
  void getMinuitValues () const;
  TMinuit* minuit; 
private:
  double overrideCallLimit; 
};

#endif 
