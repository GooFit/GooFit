#include "goofit/PdfBase.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"
#include "goofit/fitting/Params.h"

namespace GooFit {

Params::Params(PdfBase &pdf)
    : pdf_(&pdf) {
    vars_ = pdf_->getParameters();

    for(Variable *var : vars_) {
        bool added;
        if(var->IsFixed()) {
            added = Add(var->getName(), var->getValue());
        } else if(var->getLowerLimit() == var->getUpperLimit()) {
            added = Add(var->getName(), var->getValue(), var->getError());
        } else {
            added = Add(var->getName(), var->getValue(), var->getError(), var->getLowerLimit(), var->getUpperLimit());
        }

        if(!added)
            throw std::runtime_error("The name " + var->getName() + " appears more than once!");

        var->setFitterIndex(Index(var->getName()));
    }
}

void Params::SetGooFitParams(const Minuit2::MnUserParameterState &input) {
    for(Variable *var : vars_) {
        size_t counter = var->getFitterIndex();
        var->setValue(input.Value(counter));
        var->setError(input.Error(counter));
        SetValue(counter, var->getValue());
        SetError(counter, var->getError());
    }
}

} // namespace GooFit
