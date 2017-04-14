#pragma once

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameters.h>
#include <vector>
#include <memory>

class PdfBase;
class Variable;

namespace Minuit2 = ROOT::Minuit2;

namespace GooFit {

class FCN : public Minuit2::FCNBase {
protected:
    PdfBase* pdfPointer;
    std::vector<Variable*> vars;
    int numPars;
    std::unique_ptr<Minuit2::MnUserParameters> params;
    
public:
    // Access to the paramters
    Minuit2::MnUserParameters* get_params() {
        return params.get();
    }
    
    // This value is 0.5 for ll, 1 for chi2
    double Up() const override {
        return 0.5;
    }
    
    FCN(PdfBase* dat);
    
    // Minimization function
    double operator()(const std::vector<double>& pars) const override;
};

}

