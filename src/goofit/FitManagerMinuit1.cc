#include "goofit/PdfBase.h"
#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/PDFs/GooPdf.h"
#include <cstdio>
#include <cassert>
#include <limits>
#include <typeinfo>
#include <set>
#include "goofit/Variable.h"
#include "goofit/Log.h"

namespace GooFit {

Minuit1::Minuit1(PdfBase* pdfPointer) : TMinuit(max_index(pdfPointer->getParameters())+1), pdfPointer(pdfPointer), vars(pdfPointer->getParameters()) {
    int counter = 0;
    
    for(Variable* var : vars) {
        var->setFitterIndex(counter);
        DefineParameter(counter,
                                var->name.c_str(),
                                var->value,
                                var->error,
                                var->lowerlimit,
                                var->upperlimit);
        
        if(var->fixed)
            FixParameter(counter);
        
        counter++;
    }
    
    pdfPointer->copyParams();

}


Int_t Minuit1::Eval(
            int npar, 
            double* gin,
            double& fun, 
            double* fp,  
            int iflag) {
    
    std::vector<double> pars {fp, fp+npar};
    
    for(Variable* var : vars) {
        if(std::isnan(pars.at(var->getFitterIndex())))
            std::cout << "Variable " << var->name << " " << var->getIndex() << " is NaN\n";
        
        var->unchanged_ = var->value == pars.at(var->getFitterIndex());
        var->value = pars.at(var->getIndex()); //  + var->blind
    }
    
    pdfPointer->copyParams(pars);
    
    GOOFIT_TRACE("Calculating NLL");
    fun = pdfPointer->calculateNLL();
    host_callnumber++;
    return 0;
}

void FitManagerMinuit1::fit() {
    host_callnumber = 0;
    
    for(Variable* var : minuit_.getVaraibles())
        var->unchanged_ = false;

    if(0 < overrideCallLimit) {
        std::cout << "Calling MIGRAD with call limit " << overrideCallLimit << std::endl;
        double plist[1] = {overrideCallLimit};
        int err = 0;

        if(_useHesseBefore)
            minuit_.mnexcm("HESSE", plist, 1, err);

        minuit_.mnexcm("MIGRAD", plist, 1, err);

        if(_useHesse)
            minuit_.mnexcm("HESSE", plist, 1, err);

        if(_useMinos)
            minuit_.mnexcm("MINOS", plist, 1, err);

        if(_useImprove)
            minuit_.mnexcm("IMPROVE", plist, 1, err);
    } else
        minuit_.Migrad();
    
    for(Variable* var : minuit_.getVaraibles())
        minuit_.GetParameter(var->getFitterIndex(), var->value, var->error);

}


void FitManagerMinuit1::getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) {
    minuit_.mnstat(fmin, fedm, errdef, npari, nparx, istat);
    std::cout << "mnstat(fmin = " << fmin << ", fedm = " << fedm << ", errdef = " << errdef
              << ", npari = " << npari << ", nparx = " << nparx << ", istat = " << istat << ")" << std::endl;
}

}



