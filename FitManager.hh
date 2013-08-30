#ifndef PDFFUNCTOR_HH
#define PDFFUNCTOR_HH

#include "GlobalCudaDefines.hh" 
#include "EngineCore.hh"

// Glue class that talks to MINUIT
#define MINUIT_VERSION 1

#if MINUIT_VERSION == 1
#include "FitManagerMinuit1.hh"
#elif MINUIT_VERSION == 2
#include "FitManagerMinuit2.hh"
#else
#include "FitManagerMinuit3.hh"
#endif
#endif 
