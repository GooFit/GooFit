#ifndef PDFFUNCTOR_MINUIT3_HH
#define PDFFUNCTOR_MINUIT3_HH

#include "TVirtualFitter.h" 
extern PdfBase* pdfPointer; 
extern int numPars; 
#ifdef OMP_ON
#pragma omp threadprivate (numPars)
#pragma omp threadprivate (pdfPointer)
#endif
void FitFun (int &npar, double *gin, double &fun, double *fp, int iflag); 

class FitManager {
  friend class FitManager; 
public:
  FitManager (PdfBase* dat); 
  void fit (); 
  void getMinuitValues () const;

private:
  TVirtualFitter* fitter; 
};

#endif
