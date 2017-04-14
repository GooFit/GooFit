#include "goofit/PdfBase.h"
#include "goofit/FCN.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

FCN::FCN(PdfBase* dat) : pdfPointer(dat) {
    host_callnumber = 0;
    params.reset(new ROOT::Minuit2::MnUserParameters());
    pdfPointer->getParameters(vars);
    numPars = vars.size();
    
    int maxIndex = 0;

    for(Variable* var : vars) {
        if(var->lowerlimit == var->upperlimit)
            params->Add(var->name, var->value, var->error);
        else
            params->Add(var->name, var->value, var->error, var->lowerlimit, var->upperlimit);

        if(var->fixed)
            params->Fix(params->Index(var->name));

        if(maxIndex < var->getIndex())
            maxIndex = var->getIndex();
    }

    numPars = maxIndex+1;
}

double FCN::operator()(const vector<double>& pars) const {
    vector<double> gooPars; // Translates from Minuit indexing to GooFit indexing
    gooPars.resize(numPars);
    int counter = 0;

    for(Variable* var : vars) {
        gooPars[var->index] = pars[counter++];
    }

    pdfPointer->copyParams(gooPars);
    double nll = pdfPointer->calculateNLL();
    host_callnumber++;

    return nll;
}

}
