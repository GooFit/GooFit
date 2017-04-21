#include "goofit/PdfBase.h"
#include "goofit/FitManager.h"
#include "goofit/PDFs/GooPdf.h"
#include <cstdio>
#include <cassert>
#include <limits>
#include <typeinfo>
#include <set>
#include "goofit/Variable.h"

PdfBase* pdfPointer;
FitManager* currGlue = 0;
int numPars = 0;
vector<Variable*> vars;

void specialTddpPrint(double fun);

namespace GooFit {

FitManagerMinuit1::FitManagerMinuit1(PdfBase* dat)
    : minuit(0)
    , overrideCallLimit(-1)
    , _useHesseBefore(true)
    , _useHesse(true)
    , _useMinos(false)
    , _useImprove(false) {
    pdfPointer = dat;
    currGlue = this;
}

FitManagerMinuit1::~FitManagerMinuit1() {
    if(minuit)
        delete minuit;
}

void FitManagerMinuit1::setupMinuit() {
    vars.clear();
    pdfPointer->getParameters(vars);

    numPars = vars.size();

    if(minuit)
        delete minuit;

    minuit = new TMinuit(numPars);
    int maxIndex = 0;
    int counter = 0;

    for(std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
        minuit->DefineParameter(counter, (*i)->name.c_str(), (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit);

        if((*i)->fixed)
            minuit->FixParameter(counter);

        counter++;

        if(maxIndex < (*i)->getIndex())
            maxIndex = (*i)->getIndex();
    }

    numPars = maxIndex+1;
    pdfPointer->copyParams();
    minuit->SetFCN(FitFun);
}

void FitManagerMinuit1::fit() {
    setupMinuit();
    runMigrad();
}

void FitManagerMinuit1::runMigrad() {
    assert(minuit);
    host_callnumber = 0;

    if(0 < overrideCallLimit) {
        std::cout << "Calling MIGRAD with call limit " << overrideCallLimit << std::endl;
        double plist[1];
        plist[0] = overrideCallLimit;
        int err = 0;

        if(_useHesseBefore)
            minuit->mnexcm("HESSE", plist, 1, err);

        minuit->mnexcm("MIGRAD", plist, 1, err);

        if(_useHesse)
            minuit->mnexcm("HESSE", plist, 1, err);

        if(_useMinos)
            minuit->mnexcm("MINOS", plist, 1, err);

        if(_useImprove)
            minuit->mnexcm("IMPROVE", plist, 1, err);
    } else
        minuit->Migrad();
}

void FitManagerMinuit1::getMinuitValues() const {
    int counter = 0;

    for(std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
        minuit->GetParameter(counter++, (*i)->value, (*i)->error);
    }
}

void FitManagerMinuit1::getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) const {
    minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    std::cout << "mnstat(fmin = " << fmin << ", fedm = " << fedm << ", errdef = " << errdef
              << ", npari = " << npari << ", nparx = " << nparx << ", istat = " << istat << ")" << std::endl;
}

}

void FitFun(int& npar, double* gin, double& fun, double* fp, int iflag) {
    vector<double> pars;
    // Notice that npar is number of variable parameters, not total.
    pars.resize(numPars);
    int counter = 0;

    for(std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
        if(std::isnan(fp[counter]))
           std::cout << "Variable " << (*i)->name << " " << (*i)->index << " is NaN\n";

        pars[(*i)->getIndex()] = fp[counter++] + (*i)->blind;
    }

    pdfPointer->copyParams(pars);
    fun = pdfPointer->calculateNLL();
    host_callnumber++;

}


