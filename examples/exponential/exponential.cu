#include "Variable.hh" 
#include "FitManager.hh"
#include "UnbinnedDataSet.hh" 
#include "ExpPdf.hh" 
#include <iostream>

using namespace std; 

int main (int argc, char** argv) {
  // Independent variable. 
  Variable* xvar = new Variable("xvar", 0, log(RAND_MAX)); 
  
  // Data set
  UnbinnedDataSet data(xvar);
  // Generate toy events. 
  for (int i = 0; i < 100000; ++i) {
    xvar->value = xvar->upperlimit - log(1+rand());
    if (xvar->value < 0) continue;
    data.addEvent(); 
  }
  
  // Fit parameter
  Variable* alpha = new Variable("alpha", -2, 0.1, -10, 10);
  // GooPdf object
  ExpPdf* exppdf = new ExpPdf("exppdf", xvar, alpha); 
  exppdf->setData(&data);

  FitManager fitter(exppdf);
  fitter.fit(); 

  return 0;
}
