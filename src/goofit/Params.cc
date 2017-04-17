#include "goofit/PdfBase.h"
#include "goofit/fitting/Params.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

    Params::Params(PdfBase &pdf) : pdf_(&pdf) {
    pdf_->getParameters(vars_);
    num_ = vars_.size();
    
    int maxIndex = 0;
    
    for(Variable* var : vars_) {
        if(var->lowerlimit == var->upperlimit)
            Add(var->name, var->value, var->error);
        else
            Add(var->name, var->value, var->error, var->lowerlimit, var->upperlimit);
        
        if(var->fixed)
            Fix(Index(var->name));
        
        if(maxIndex < var->getIndex())
            maxIndex = var->getIndex();
    }
    
    num_ = maxIndex+1;
}
    
void Params::GetParams() const {
    int counter = 0;
    for(Variable* var : vars_) {
        var->value = Value(counter);
        var->error = Error(counter++);
    }
}
    
}
