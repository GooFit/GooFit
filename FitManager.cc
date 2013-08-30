#include "PdfBase.hh"
#include "FitManager.hh"
#include "GooPdf.hh" 
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
#include "FitManagerMinuit1.cc"
#elif MINUIT_VERSION == 2
#include "FitManagerMinuit2.cc"
#else 
#include "FitManagerMinuit3.cc"
#endif 
