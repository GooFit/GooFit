#include "goofit/PdfBase.h"
#include "goofit/fitting/FCN.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

FCN::FCN(Params& params) : params_(&params) {
    host_callnumber = 0;
    
}

double FCN::operator()(const vector<double>& pars) const {
    vector<double> gooPars; // Translates from Minuit indexing to GooFit indexing
    gooPars.resize(params_->num_);
    int counter = 0;

    for(Variable* var : params_->vars_) {
        gooPars[var->index] = pars[counter++];
    }

    params_->pdf_->copyParams(gooPars);
    double nll = params_->pdf_->calculateNLL();
    host_callnumber++;

    return nll;
}

}
