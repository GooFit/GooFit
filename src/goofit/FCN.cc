#include "goofit/PdfBase.h"
#include "goofit/Log.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"
#include "goofit/fitting/FCN.h"

namespace GooFit {

FCN::FCN(Params &params)
    : params_(&params) {
    host_callnumber = 0;

    // Verify that all varaibles need to be recached
    for(Variable *var : params_->vars_)
        var->setChanged(true);
}

double FCN::operator()(const std::vector<double> &pars) const {
    // Translate from Minuit indexing to GooFit indexing
    std::vector<double> gooPars(max_index(params_->vars_) + 1);

    for(Variable *var : params_->vars_) {
        var->setChanged(var->getValue() != pars.at(var->getFitterIndex()));
        gooPars.at(var->getIndex()) = pars.at(var->getFitterIndex()) - var->blind;
    }

    params_->pdf_->copyParams(gooPars);

    GOOFIT_TRACE("Calculating NLL");
    double nll = params_->pdf_->calculateNLL();

    host_callnumber++;

    return nll;
}

double FCN::operator()() const {
    std::vector<double> pars = makePars();

    // Translate from Minuit indexing to GooFit indexing
    std::vector<double> gooPars(max_index(params_->vars_) + 1);

    for(Variable *var : params_->vars_) {
        var->setChanged(true);
        gooPars.at(var->getIndex()) = pars.at(var->getFitterIndex()) - var->blind;
    }

    params_->pdf_->copyParams(gooPars);

    GOOFIT_TRACE("Calculating NLL");
    double nll = params_->pdf_->calculateNLL();

    host_callnumber++;

    return nll;
}

std::vector<double> FCN::makePars() const {
    std::vector<double> minuitPars(max_fitter_index(params_->vars_) + 1);
    for(Variable *var : params_->vars_) {
        minuitPars.at(var->getFitterIndex()) = var->getValue();
    }
    return minuitPars;
}

// Get the number of variable parameters
Params *FCN::GetParams() { return params_; }

} // namespace GooFit
