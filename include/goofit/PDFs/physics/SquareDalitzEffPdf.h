//#ifndef SQUAREDALITZEFF_PDF_HH
//#define SQUAREDALITZEFF_PDF_HH

#include <goofit/PDFs/GooPdf.h>

class SquareDalitzEffPdf : public GooPdf {

public:
  // Very specific efficiency parametrisation for semileptonically-tagged D0->KSPiPi decays as determined from data
  // Uses variables of square Dalitz plot - m' and theta' 
  SquareDalitzEffPdf (std::string n, vector<Variable*> obses, vector<Variable*> coeffs, vector<Variable*> constvals); 

private:

};

#endif
