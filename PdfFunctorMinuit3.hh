#ifndef PDFFUNCTOR_MINUIT3_HH
#define PDFFUNCTOR_MINUIT3_HH

#include "TVirtualFitter.h" 
extern FunctorBase* pdfPointer; 
extern int numPars; 
#ifdef OMP_ON
#pragma omp threadprivate (numPars)
#pragma omp threadprivate (pdfPointer)
#endif
void FitFun (int &npar, double *gin, double &fun, double *fp, int iflag); 

class PdfFunctor {
  friend class PdfFunctor; 
public:
  PdfFunctor (FunctorBase* dat); 
  void fit (); 
  void getMinuitValues () const;

private:
  TVirtualFitter* fitter; 
};

#endif
