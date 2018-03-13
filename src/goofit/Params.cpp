#include <goofit/PDFs/GooPdf.h>
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>
#include <goofit/fitting/Params.h>

namespace GooFit {

Params::Params(PdfBase &pdf)
    : pdf_(&pdf) {
    vars_ = pdf_->getParameters();

    for(Variable &var : vars_) {
        bool added;
        if(var.IsFixed()) {
            added = Add(var.getName(), var.getValue());
        } else if(var.getLowerLimit() == var.getUpperLimit()) {
            added = Add(var.getName(), var.getValue(), var.getError());
        } else {
            added = Add(var.getName(), var.getValue(), var.getError(), var.getLowerLimit(), var.getUpperLimit());
        }

        if(!added)
            throw std::runtime_error("The name " + var.getName() + " appears more than once!");

        var.setFitterIndex(Index(var.getName()));
    }
}

void Params::SetGooFitParams(const Minuit2::MnUserParameterState &input) {
    for(Variable &var : vars_) {
        int counter = var.getFitterIndex();
        var.setValue(input.Value(counter));
        var.setError(input.Error(counter));
        SetValue(counter, var.getValue());
        SetError(counter, var.getError());
    }
}

std::vector<double> Params::make_minuit_vector() const {
    std::vector<double> minuitPars(max_fitter_index(vars_) + 1);
    for(const Variable &var : vars_) {
        minuitPars.at(var.getFitterIndex()) = var.getValue();
    }
    return minuitPars;
}

void Params::from_minuit_vector(const std::vector<double> &values, bool force_changed) {
    for(Variable &var : vars_) {
        var.setChanged(force_changed ? true : var.getValue() != values.at(var.getFitterIndex()));
        pdf_->updateVariable(var, values.at(var.getFitterIndex()) - var.getBlind(Variable::Key()));
    }

    if(do_record_)
        recorded_.push_back(values);

    pdf_->updateParameters();
}

} // namespace GooFit
