#include "goofit/PdfBase.h"
#include "goofit/fitting/Params.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

    Params::Params(PdfBase &pdf) : pdf_(&pdf) {
    vars_ = pdf_->getParameters();
        
    for(Variable* var : vars_) {
        bool added;
        if(var->IsFixed()) {
            added = Add(var->GetName(), var->GetValue());
        } else if(var->GetLowerLimit() == var->GetUpperLimit()) {
            added = Add(var->GetName(), var->GetValue(), var->GetError());
        } else {
            added = Add(var->GetName(), var->GetValue(), var->GetError(), var->GetLowerLimit(), var->GetUpperLimit());
        }
        
        if(!added)
            throw std::runtime_error("The name " + var->GetName() + " appears more than once!");
        
        var->SetFitterIndex(Index(var->GetName()));
    }
    
}
    
void Params::SetGooFitParams(const Minuit2::MnUserParameterState& input) {
    for(Variable* var : vars_) {
        size_t counter = var->GetFitterIndex();
        var->SetValue(input.Value(counter));
        var->SetError(input.Error(counter));
        SetValue(counter, var->value);
        SetError(counter, var->error);
    }
}
    
}
