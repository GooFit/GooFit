#include "goofit/PdfBase.h"
#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/PDFs/GooPdf.h"
#include <cstdio>
#include <cassert>
#include <limits>
#include <typeinfo>
#include <set>
#include "goofit/Variable.h"

PdfBase* pdfPointer;
std::vector<Variable*> vars;

namespace GooFit {

FitManagerMinuit1::FitManagerMinuit1(PdfBase* ptr) {
    pdfPointer = ptr;
    pdfPointer->getParameters(vars);

    minuit_.reset(new TMinuit(max_index(vars)+1));

    int counter = 0;

    for(Variable* var : vars) {
        var->setFitterIndex(counter);
        minuit_->DefineParameter(counter,
                                 var->name.c_str(),
                                 var->value,
                                 var->error,
                                 var->lowerlimit,
                                 var->upperlimit);

        if(var->fixed)
            minuit_->FixParameter(counter);

        counter++;
    }

    pdfPointer->copyParams();
    minuit_->SetFCN(FitFun);
}

void FitManagerMinuit1::runFit() {
    assert(minuit_);
    host_callnumber = 0;
    
    for(Variable* var : vars)
        var->unchanged_ = false;

    if(0 < overrideCallLimit) {
        std::cout << "Calling MIGRAD with call limit " << overrideCallLimit << std::endl;
        double plist[1] = {overrideCallLimit};
        int err = 0;

        if(_useHesseBefore)
            minuit_->mnexcm("HESSE", plist, 1, err);

        minuit_->mnexcm("MIGRAD", plist, 1, err);

        if(_useHesse)
            minuit_->mnexcm("HESSE", plist, 1, err);

        if(_useMinos)
            minuit_->mnexcm("MINOS", plist, 1, err);

        if(_useImprove)
            minuit_->mnexcm("IMPROVE", plist, 1, err);
    } else
        minuit_->Migrad();
    
    for(Variable* var : vars)
        minuit_->GetParameter(var->getFitterIndex(), var->value, var->error);

}


void FitManagerMinuit1::getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) const {
    minuit_->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    std::cout << "mnstat(fmin = " << fmin << ", fedm = " << fedm << ", errdef = " << errdef
              << ", npari = " << npari << ", nparx = " << nparx << ", istat = " << istat << ")" << std::endl;
}

}

/// Fit function for Minuit
void FitFun(int& npar, //< THe number of parameters
            double* gin, //< The derivatives can be stored here if flag is 2 (output)
            double& fun, //< The value of the function at this point (output)
            double* fp, //< The input parameters
            int iflag //< This is 1 the first time, 2, for derivatives, and 3 after the fit is finished. It is something else if computing.
            ) {

    std::vector<double> pars {fp, fp+npar};

    for(Variable* var : vars) {
        if(std::isnan(pars.at(var->getFitterIndex())))
            std::cout << "Variable " << var->name << " " << var->getIndex() << " is NaN\n";
        
        var->unchanged_ = var->value == pars.at(var->getFitterIndex());
        pars.at(var->getIndex()) = var->value; //  + var->blind
    }

    pdfPointer->copyParams(pars);
    fun = pdfPointer->calculateNLL();
    host_callnumber++;
}

