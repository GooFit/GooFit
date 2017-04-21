
#include <cstdio>
#include <cassert>
#include <limits>
#include <typeinfo>
#include <set>

#include "goofit/PdfBase.h"
#include "goofit/FitManagerMinuit2.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

ROOT::Minuit2::FunctionMinimum* FitManagerMinuit2::fit() {
    host_callnumber = 0;
    params = new ROOT::Minuit2::MnUserParameters();
    vars.clear();
    pdfPointer->getParameters(vars);

    numPars = vars.size();
    int maxIndex = 0;

    for(std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
        if((*i)->lowerlimit == (*i)->upperlimit)
            params->Add((*i)->name, (*i)->value, (*i)->error);
        else
            params->Add((*i)->name, (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit);

        if((*i)->fixed)
            params->Fix(params->Index((*i)->name));

        if(maxIndex < (*i)->getIndex())
            maxIndex = (*i)->getIndex();
    }

    numPars = maxIndex+1;
    migrad = new ROOT::Minuit2::MnMigrad(*this, *params);
    ROOT::Minuit2::FunctionMinimum* ret = new ROOT::Minuit2::FunctionMinimum((*migrad)());

    return ret;
}

double FitManagerMinuit2::operator()(const vector<double>& pars) const {
    vector<double> gooPars; // Translates from Minuit indexing to GooFit indexing
    gooPars.resize(numPars);
    int counter = 0;

    for(std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
        gooPars[(*i)->index] = pars[counter++];
    }

    pdfPointer->copyParams(gooPars);
    double nll = pdfPointer->calculateNLL();
    host_callnumber++;

    return nll;
}
}
