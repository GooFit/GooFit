#pragma once

#include <memory>
#include <TMinuit.h>
#include <Minuit2/FunctionMinimum.h>
#include "goofit/fitting/FCN.h"

class PdfBase;

namespace GooFit {

class FitManagerMinuit2 {
public:
    FitManagerMinuit2(PdfBase* dat);
    
    // This runs the fit
    ROOT::Minuit2::FunctionMinimum fit();
    
private:
    Params upar_;
    FCN fcn_;
};
    
}
