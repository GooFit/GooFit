#include "goofit/PdfBase.h"
#include "goofit/fitting/Params.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

#include <unordered_set>

namespace GooFit {

    Params::Params(PdfBase &pdf) : pdf_(&pdf) {
    pdf_->getParameters(vars_);
        
    for(Variable* var : vars_) {
        bool added;
        if(var->fixed) {
            added = Add(var->name, var->value);
        } else if(var->lowerlimit == var->upperlimit) {
            added = Add(var->name, var->value, var->error);
        } else {
            added = Add(var->name, var->value, var->error, var->lowerlimit, var->upperlimit);
        }
        
        if(!added)
            throw std::runtime_error("The name " + var->name + " appears more than once!");
        
        var->index = Index(var->name);
    }
    
}
    
void Params::SetGooFitParams(const Minuit2::MnUserParameterState& input) {
    for(Variable* var : vars_) {
        size_t counter = Index(var->name);
        var->value = input.Value(counter);
        var->error = input.Error(counter);
        SetValue(counter, var->value);
        SetError(counter, var->error);
    }
}
    
}
