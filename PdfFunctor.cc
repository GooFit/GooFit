#include "FunctorBase.hh"
#include "PdfFunctor.hh"
#include "EngineCore.hh" 
#include <cstdio> 
#include <cassert> 
#include <limits> 
#include <typeinfo> 
#include <set>
#include "Variable.hh" 

using namespace std; 

#ifdef OMP_ON
#pragma omp threadprivate (host_params)
#endif

#if MINUIT_VERSION == 1
#include "PdfFunctorMinuit1.cc"
#elif MINUIT_VERSION == 2
#include "PdfFunctorMinuit2.cc"
#else 
#include "PdfFunctorMinuit3.cc"
#endif 
