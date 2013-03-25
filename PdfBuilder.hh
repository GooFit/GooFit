#ifndef PDFBUILDER_HH
#define PDFBUILDER_HH

#include "GlobalCudaDefines.hh" 
#include "ThrustPdfFunctor.hh"

// Glue class that talks to MINUIT

#define MINUIT_VERSION 1

#if MINUIT_VERSION == 2
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
class PdfFunctor : public ROOT::Minuit2::FCNBase {
  friend class PdfBuilder; 
public:
  PdfFunctor (FunctorBase* dat) {pdfPointer = dat;} 
  virtual double Up() const {return 1.0;}
  double operator () (const std::vector<double>& pars) const; 
  ROOT::Minuit2::FunctionMinimum* fit (); 
  void getMinuitValues () const;

protected:
  FunctorBase* pdfPointer; 
  ROOT::Minuit2::MnUserParameters* params;
  ROOT::Minuit2::MnMigrad* migrad; 

  std::set<Variable*> vars; 
  int numPars; 
};
#elif MINUIT_VERSION == 3
#include "TVirtualFitter.h" 
extern FunctorBase* pdfPointer; 
extern int numPars; 
#ifdef OMP_ON
#pragma omp threadprivate (numPars)
#pragma omp threadprivate (pdfPointer)
#endif
void FitFun (int &npar, double *gin, double &fun, double *fp, int iflag); 

class PdfFunctor {
  friend class PdfBuilder; 
public:
  PdfFunctor (FunctorBase* dat); 
  void fit (); 
  void getMinuitValues () const;

private:
  TVirtualFitter* fitter; 
};
#else

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
  void fit (); 
  void getMinuitValues () const;
private:
  TMinuit* minuit; 
  double overrideCallLimit; 
};

#endif 

#endif
