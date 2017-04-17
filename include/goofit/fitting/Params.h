#pragma once

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>
#include <vector>

class PdfBase;
class Variable;

namespace Minuit2 = ROOT::Minuit2;

namespace GooFit {

class FCN;
    
class Params : public Minuit2::MnUserParameters {
    friend FCN;
protected:
    std::vector<Variable*> vars_;
    PdfBase* pdf_;
    size_t num_;
    
public:
    using MnUserParameters::MnUserParameters;
    
    Params(PdfBase &pdf);
    
    // Read the values back into GooFit
    void SetGooFitParams(const Minuit2::MnUserParameterState& input);
    
    // Get the number of params in the fit
    size_t GetNumPars() const {return num_;};
};
    
}

