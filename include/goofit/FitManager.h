#ifndef FITMANAGER_HH
#define FITMANAGER_HH

#include "goofit/GlobalCudaDefines.h"
#include "goofit/PDFs/GooPdf.h"

// Glue class that talks to MINUIT
#define MINUIT_VERSION 1

#if MINUIT_VERSION == 1
#include "goofit/FitManagerMinuit1.h"
#elif MINUIT_VERSION == 2
#include "goofit/FitManagerMinuit2.h"
#else
#include "goofit/FitManagerMinuit3.h"
#endif
#endif
