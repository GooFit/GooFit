#pragma once

#include <memory>
#include <TMinuit.h>
#include <Minuit2/FunctionMinimum.h>
#include "goofit/FCN.h"

class PdfBase;
namespace GooFit {

class FitManagerMinuit2 {
public:
    FitManagerMinuit2(PdfBase* dat);
    
    // This runs the fit
    ROOT::Minuit2::FunctionMinimum fit();
    
    
private:
    std::unique_ptr<FCN> fcn_;
};
    
}
