#include "goofit/fitting/FitManagerMinuit2.h"

#include "goofit/Color.h"

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include "Minuit2/MnPrint.h"

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
    std::cout << GooFit::gray << GooFit::bold;
    Minuit2::FunctionMinimum min = migrad(maxfcn_);
    
    // Print nice output
    std::cout << GooFit::reset << (min.IsValid() ? GooFit::green : GooFit::red);   
    std::cout << min << GooFit::reset;
    std::cout << GooFit::magenta << timer << GooFit::reset << std::endl;
    
    if(min.IsValid()) {
        retval_ = FitErrors::Valid;
    } else {
        std::cout << GooFit::red;
        std::cout << "HesseFailed: " << min.HesseFailed() << std::endl;
        std::cout << "HasCovariance: " << min.HasCovariance() << std::endl;
        std::cout << "HasValidCovariance: " << min.HasValidCovariance() << std::endl;
        std::cout << "HasValidParameters: " << min.HasValidParameters() << std::endl;
        std::cout << "IsAboveMaxEdm: " << min.IsAboveMaxEdm() << std::endl;
        std::cout << "HasReachedCallLimit: " << min.HasReachedCallLimit() << std::endl;
        std::cout << "HasAccurateCovar: " << min.HasAccurateCovar() << std::endl;
        std::cout << "HasPosDefCovar : " << min.HasPosDefCovar () << std::endl;
        std::cout << "HasMadePosDefCovar : " << min.HasMadePosDefCovar () << std::endl;
        std::cout << GooFit::reset;
        
        retval_ = FitErrors::InValid;
    }
    
    // Set the parameters in GooFit to the new values
    upar_.SetGooFitParams(min.UserState());

    Minuit2::MnPrint::SetLevel(val);
    return min;
}
    
}
  

