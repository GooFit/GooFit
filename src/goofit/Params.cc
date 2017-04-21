#include "goofit/PdfBase.h"
#include "goofit/fitting/Params.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

#include <unordered_set>

namespace GooFit {

    Params::Params(PdfBase &pdf) : pdf_(&pdf) {
    pdf_->getParameters(vars_);
    num_ = vars_.size();
    
    size_t maxIndex = 0;
        
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
        
        if(maxIndex < var->getIndex())
            maxIndex = var->getIndex();
    }
    
    num_ = maxIndex+1;
}
    
void Params::SetGooFitParams(const Minuit2::MnUserParameterState& input) {
    int counter = 0;
    for(Variable* var : vars_) {
        var->value = input.Value(counter);
        var->error = input.Error(counter);
        SetValue(counter, var->value);
        SetError(counter, var->error);
        counter++;
    }
}
    
}
