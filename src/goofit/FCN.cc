#include "goofit/PdfBase.h"
#include "goofit/fitting/FCN.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"
#include "goofit/Log.h"

namespace GooFit {

FCN::FCN(Params& params) : params_(&params) {
    host_callnumber = 0;
    
    // Verify that all varaibles need to be recached
    for(Variable* var : params_->vars_)
        var->SetChanged(true);
    
}

double FCN::operator()(const std::vector<double>& pars) const {
    
    // Translate from Minuit indexing to GooFit indexing
    std::vector<double> gooPars;
    gooPars.resize(max_index(params_->vars_)+1);
    
    for(Variable* var : params_->vars_) {
        // TODO: support for blinding
        var->SetChanged(var->value != pars.at(var->GetFitterIndex()));
        gooPars.at(var->GetIndex()) = pars.at(var->GetFitterIndex());
    }

    params_->pdf_->copyParams(gooPars);
    
    GOOFIT_TRACE("Calculating NLL");
    double nll = params_->pdf_->calculateNLL();
    
    host_callnumber++;

    return nll;
}

// Get the number of variable parameters
Params* FCN::GetParams() {
    return params_;
}
    
}
