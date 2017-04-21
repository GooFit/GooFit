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

#ifdef PRINTCALLS
    specialTddpPrint(fun);
#endif
}


#ifdef PRINTCALLS
void specialTddpPrint(double fun) {
    // Stupid amplitude-fit debugging method.
   std::cout << "Function call " << host_callnumber << ": " << fun << "\n";
    currGlue->getMinuitValues();
    int varCount = 1;

    for(std::vector<Variable*>::iterator v = vars.begin(); v != vars.end(); ++v) {
        if(!(*v))
           std::cout << "Null!" <<std::endl;

        if((*v)->fixed)
            continue;

        const fptype _mD0 = 1.86484;
        const fptype _mD02 = _mD0 *_mD0;
        const fptype _mD02inv = 1./_mD02;
        double stupidSpecialModifier = 1; // Mikhail interprets some of the weights differently.

        if(((*v)->name == "f0_980_amp_real") ||
                ((*v)->name == "f0_980_amp_imag") ||
                ((*v)->name == "f0_1370_amp_real") ||
                ((*v)->name == "f0_1370_amp_imag") ||
                ((*v)->name == "f0_1500_amp_real") ||
                ((*v)->name == "f0_1500_amp_imag") ||
                ((*v)->name == "f0_1710_amp_real") ||
                ((*v)->name == "f0_1710_amp_imag") ||
                ((*v)->name == "f0_600_amp_real") ||
                ((*v)->name == "f0_600_amp_imag"))
            stupidSpecialModifier = -_mD02;
        else if(((*v)->name == "f2_1270_amp_real") ||
                ((*v)->name == "f2_1270_amp_imag"))
            stupidSpecialModifier = -_mD02inv;
        else if(((*v)->name == "nonr_amp_real") ||
                ((*v)->name == "nonr_amp_imag"))
            stupidSpecialModifier = -1;

       std::cout.width(20);
       std::cout << (*v)->name;
       std::cout.setf(ios_base::right, ios_base::adjustfield);
       std::cout.width(3);
       std::cout << varCount++;
       std::cout.setf(ios_base::right, ios_base::adjustfield);
       std::cout.precision(8);
       std::cout << "  ";
       std::cout.width(12);
       std::cout << (*v)->value / stupidSpecialModifier;
       std::cout.setf(ios_base::right, ios_base::adjustfield);
       std::cout.precision(8);
       std::cout << "  ";
       std::cout.width(12);
       std::cout << (*v)->error;
       std::cout <<std::endl;
    }

   std::cout <<std::endl;
}
#endif

