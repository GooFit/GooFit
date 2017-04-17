#include "goofit/PdfBase.h"
#include "goofit/fitting/FCN.h"
#include "goofit/fitting/FitManagerMinuit2.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/FunctionMinimum.h>
#include "Minuit2/MnPrint.h"

#include <rang.hpp>
#include <CLI/Timer.hpp>

namespace GooFit {
    
FitManagerMinuit2::FitManagerMinuit2(PdfBase* dat) : upar_(*dat), fcn_(upar_) {}

Minuit2::FunctionMinimum FitManagerMinuit2::fit() {
    auto val = Minuit2::MnPrint::Level();
    Minuit2::MnPrint::SetLevel(3);
    
    // Setting global call number to 0
    host_callnumber = 0;
    
    CLI::Timer timer{"The minimization took"};
    
    Minuit2::MnMigrad migrad{fcn_, upar_};
    
    // Do the minimization
    std::cout << rang::fg::gray << rang::style::bold;
    Minuit2::FunctionMinimum min = migrad(maxfcn_);
    
    // Print nice output
    std::cout << rang::style::reset << (min.IsValid() ? rang::fg::green : rang::fg::red);   
    std::cout << min << rang::style::reset;
    std::cout << rang::fg::magenta << timer << rang::style::reset << std::endl;
    
    // Set the parameters in GooFit to the new values
    upar_.SetGooFitParams(min.UserState());

    Minuit2::MnPrint::SetLevel(val);
    return min;
}
    
}
  

