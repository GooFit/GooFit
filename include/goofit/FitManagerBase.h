#pragma once

#include "goofit/GlobalCudaDefines.h"
#include <vector>

class PdfBase;
class Variable;

/// Base class for FitManagers
class FitManagerBase {

public:
    virtual ~FitManagerBase() = default;
    
    // This is expected to prepare the fitting
    virtual void setup() {};
    
    // This is required to run the fit
    virtual void runFit() = 0;
    
    // This runs the fit
    virtual void fit() {
        setup();
        runFit();
    };
    
};
