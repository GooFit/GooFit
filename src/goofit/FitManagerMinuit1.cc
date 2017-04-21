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
vector<Variable*> vars;
int numPars = 0;

namespace GooFit {

void FitManagerMinuit1::setup() {
    vars.clear();
    pdfPointer->getParameters(vars);

    // Total number of params
    numPars = vars.size();

    minuit_.reset(new TMinuit(numPars));

    int maxIndex = 0;
    int counter = 0;

    for(Variable* var : vars) {
        minuit_->DefineParameter(counter,
                                 var->name.c_str(),
                                 var->value,
                                 var->error,
                                 var->lowerlimit,
                                 var->upperlimit);

        if(var->fixed)
            minuit_->FixParameter(counter);

        counter++;

        if(maxIndex < var->getIndex())
            maxIndex = var->getIndex();
    }

    numPars = maxIndex+1;
    pdfPointer->copyParams();
    minuit_->SetFCN(FitFun);
}

void FitManagerMinuit1::runFit() {
    assert(minuit_);
    host_callnumber = 0;

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
}

void FitManagerMinuit1::getMinuitValues() const {
    int counter = 0;

    for(Variable* var : vars) {
        minuit_->GetParameter(counter++, var->value, var->error);
    }
}

void FitManagerMinuit1::getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) const {
    minuit_->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    std::cout << "mnstat(fmin = " << fmin << ", fedm = " << fedm << ", errdef = " << errdef
              << ", npari = " << npari << ", nparx = " << nparx << ", istat = " << istat << ")" << std::endl;
}

}

void FitFun(int& npar, double* gin, double& fun, double* fp, int iflag) {
    std::vector<double> pars;
    // Notice that npar is number of variable parameters, not total.
    pars.resize(numPars);
    int counter = 0;

    for(Variable* var : vars) {
        if(std::isnan(fp[counter]))
            cout << "Variable " << var->name << " " << var->index << " is NaN\n";

        pars.at(var->getIndex()) = fp[counter++] + var->blind;
    }

    pdfPointer->copyParams(pars);
    fun = pdfPointer->calculateNLL();
    host_callnumber++;
}



