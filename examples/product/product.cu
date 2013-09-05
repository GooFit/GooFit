#include "Variable.hh" 
#include "FitManager.hh"
#include "UnbinnedDataSet.hh" 
#include "ExpPdf.hh" 
#include "ProdPdf.hh"
#include <iostream>

using namespace std; 

int main (int argc, char** argv) {
  Variable* xvar = new Variable("xvar", 0, log(1+RAND_MAX/2)); 
  Variable* yvar = new Variable("yvar", 0, log(1+RAND_MAX/2)); 
  
  vector<Variable*> varList;
  varList.push_back(xvar);
  varList.push_back(yvar);
  UnbinnedDataSet data(varList);
  for (int i = 0; i < 100000; ++i) {
    xvar->value = xvar->upperlimit - log(1+rand()/2);
    yvar->value = yvar->upperlimit - log(1+rand()/2);
    data.addEvent(); 
  }
  
  Variable* alpha_x = new Variable("alpha_x", -2.4, 0.1, -10, 10);
  Variable* alpha_y = new Variable("alpha_y", -1.1, 0.1, -10, 10);
  vector<PdfBase*> pdfList;
  pdfList.push_back(new ExpPdf("exp_x", xvar, alpha_x)); 
  pdfList.push_back(new ExpPdf("exp_y", yvar, alpha_y)); 

  ProdPdf* product = new ProdPdf("product", pdfList);
  product->setData(&data);

  FitManager fitter(product);
  fitter.fit();

  return 0;
}
