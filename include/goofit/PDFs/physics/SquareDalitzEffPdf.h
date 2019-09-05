
#pragma once


#include <vector>

#include <goofit/GlobalCudaDefines.h> // Need this for 'fptype'
#include <goofit/Variable.h>
#include <goofit/detail/Complex.h>
#include <goofit/PDFs/GooPdf.h>


namespace GooFit{

class SquareDalitzEffPdf : public GooPdf {

public:
  // Very specific efficiency parametrisation for semileptonically-tagged D0->KSPiPi decays as determined from data
  // Uses variables of square Dalitz plot - m' and theta' 
    SquareDalitzEffPdf (std::string n, std::vector<Variable*> obses, std::vector<Variable*> coeffs, std::vector<Variable*> constvals); 

private:

};

}

