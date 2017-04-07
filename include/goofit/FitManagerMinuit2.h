#pragma once

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnContours.h>
#include <Minuit2/MnPlot.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/ContoursError.h>

namespace GooFit {

class FitManagerMinuit2 : public ROOT::Minuit2::FCNBase {
    // Not needed? friend class FitManager;
public:
    FitManagerMinuit2(PdfBase* dat) {
        pdfPointer = dat;
    }
    virtual double Up() const {
        return 1.0;
    }
    double operator()(const std::vector<double>& pars) const;
    ROOT::Minuit2::FunctionMinimum* fit();
    void getMinuitValues() const;

protected:
    PdfBase* pdfPointer;
    ROOT::Minuit2::MnUserParameters* params;
    ROOT::Minuit2::MnMigrad* migrad;

    std::vector<Variable*> vars;
    int numPars;
};

}

