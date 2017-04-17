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
    
    // Set the maximum number of calls. 0 for Minuit2 default.
    void setMaxCalls(unsigned int max_calls = 0) {
        maxfcn_ = max_calls;
    }
 
    // Get a pointer to the params
    Params* getParams() {return &upar_;}
    
    // Get a pointer to the fcn
    FCN* getFCN() {return &fcn_;}
    

    // Print Minos error analysis after fitting
    void useMinos(bool use=true) {
        useMinos_ = use;
    }

    
private:
    Params upar_;
    FCN fcn_;
    unsigned int maxfcn_ {0};
    bool useMinos_ {false};
};
    
}
