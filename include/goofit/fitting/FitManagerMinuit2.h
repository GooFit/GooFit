#pragma once

#include <memory>
#include <Minuit2/FunctionMinimum.h>

#include "goofit/GlobalCudaDefines.h"
#include "goofit/fitting/FCN.h"

class PdfBase;

namespace GooFit {
    
enum class FitErrors : int {
    Valid = 0,
    NotRun = 50,
    InValid
};

class FitManagerMinuit2 {
    
public:
    FitManagerMinuit2(PdfBase* dat);
    
    /// This runs the fit
    ROOT::Minuit2::FunctionMinimum fit();
    
    /// Set the maximum number of calls. 0 for Minuit2 default.
    void setMaxCalls(unsigned int max_calls = 0) {
        maxfcn_ = max_calls;
    }
 
    /// Get a pointer to the params
    Params* getParams() {return &upar_;}
    
    /// Get a pointer to the fcn
    FCN* getFCN() {return &fcn_;}
    
    /// Check to see if fit is valid
    operator bool() {return retval_ == FitErrors::Valid;}
    
    /// Return value for program
    operator int() {return static_cast<int>(retval_);}

    
private:
    Params upar_;
    FCN fcn_;
    unsigned int maxfcn_ {0};
    FitErrors retval_ {FitErrors::NotRun};
};
    
}
