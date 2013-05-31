#ifndef PDFFUNCTOR_HH
#define PDFFUNCTOR_HH

#include "GlobalCudaDefines.hh" 
#include "ThrustPdfFunctor.hh"

// Glue class that talks to MINUIT
#define MINUIT_VERSION 1

#if MINUIT_VERSION == 1
#include "PdfFunctorMinuit1.hh"
#elif MINUIT_VERSION == 2
#include "PdfFunctorMinuit2.hh"
#else
#include "PdfFunctorMinuit3.hh"
#endif
#endif 
