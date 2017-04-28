#include "goofit/PdfBase.h"
#include "goofit/fitting/FCN.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

FCN::FCN(Params& params) : params_(&params) {
    host_callnumber = 0;
    
}

double FCN::operator()(const std::vector<double>& pars) const {
    
    // Get the max value of GooFit index
    Variable* max_ind_ptr = *std::max_element(std::begin(params_->vars_),
                                              std::end(params_->vars_),
                                              [](const Variable *a, const Variable *b)
                                                {return a->getIndex() < b->getIndex();});
    int max_ind = max_ind_ptr->getIndex();
    
    // Translate from Minuit indexing to GooFit indexing
    std::vector<double> gooPars;
    gooPars.resize(max_ind+1);
    
    for(Variable* var : params_->vars_) {
        var->unchanged_ = var->value == pars.at(var->getFitterIndex());
        gooPars.at(var->getIndex()) = pars.at(var->getFitterIndex());
    }

    params_->pdf_->copyParams(gooPars);
    double nll = params_->pdf_->calculateNLL();
    host_callnumber++;

    return nll;
}

// Get the number of variable parameters
Params* FCN::GetParams() {
    return params_;
}
    
}
