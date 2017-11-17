#include <goofit/Log.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>
#include <goofit/fitting/FCN.h>

namespace GooFit {

FCN::FCN(Params &params)
    : params_(&params) {
    host_callnumber = 0;

    // Verify that all varaibles need to be recached
    for(Variable &var : params_->vars_)
        var.setChanged(true);
}

double FCN::operator()(const std::vector<double> &pars) const {
    // Translate from Minuit indexing to GooFit indexing
    params_->from_minuit_vector(pars);

    GOOFIT_TRACE("Calculating NLL");
    double nll = params_->pdf_->calculateNLL();

    host_callnumber++;

    return nll;
}

double FCN::operator()() const {
    // Make a vector of current values, then force the recalculation step
    std::vector<double> pars = params_->make_minuit_vector();
    params_->from_minuit_vector(pars, true);

    GOOFIT_TRACE("Calculating NLL");
    double nll = params_->pdf_->calculateNLL();

    host_callnumber++;

    return nll;
}

// Get the number of variable parameters
Params *FCN::GetParams() { return params_; }

} // namespace GooFit
