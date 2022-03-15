#include <cstdio>
#include <goofit/Log.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>
#include <goofit/fitting/FitManagerMinuit1.h>
#include <limits>
#include <set>
#include <typeinfo>

namespace GooFit {

Minuit1::Minuit1(PdfBase *pdfPointer)
    : TMinuit(pdfPointer->getParameters().size() + 1)
    , pdfPointer(pdfPointer)
    , vars(pdfPointer->getParameters()) {
    int counter = 0;

    for(Variable &var : vars) {
        var.setFitterIndex(counter);

        Int_t err = DefineParameter(
            counter, var.getName().c_str(), var.getValue(), var.getError(), var.getLowerLimit(), var.getUpperLimit());

        if(GetNumPars() != counter + 1)
            throw GooFit::GeneralError(
                "Error when implementing param {} (possibly invalid error/lowerlimit/upperlimit values)!",
                var.getName());

        if(err != 0)
            throw GooFit::GeneralError("Was not able to implement param {} (error {})", var.getName(), err);

        if(var.IsFixed())
            FixParameter(counter);

        counter++;
    }

    pdfPointer->copyParams();
}

auto Minuit1::Eval(int npar, double *gin, double &fun, double *fp, int iflag) -> Int_t {
    std::vector<double> pars{fp, fp + GetNumPars()};

    std::vector<double> gooPars;
    gooPars.resize(vars.size() + 1);

    for(Variable &var : vars) {
        if(std::isnan(pars.at(var.getFitterIndex())))
            GOOFIT_WARN("Variable {} at {} is NaN", var.getName(), var.getIndex());

        var.setChanged(var.getValue() != pars.at(var.getFitterIndex()));
        var.setValue(pars.at(var.getFitterIndex()));
        // TODO: fix this to reflect Minuit2 updateVariable. Actually not needed...
        // gooPars.at(var->getIndex()) = var->getValue() - var->blind;
        pdfPointer->updateVariable(var, var.getValue() - var.getBlind(Variable::Key()));
    }

    pdfPointer->updateParameters();

    GOOFIT_TRACE("Calculating NLL");
    fun = pdfPointer->calculateNLL();
    host_callnumber++;
    return 0;
}

void FitManagerMinuit1::fit() {
    host_callnumber = 0;

    for(Variable &var : minuit_.getVaraibles())
        var.setChanged(true);

    std::cout << GooFit::gray << GooFit::bold;

    if(0 < overrideCallLimit) {
        std::cout << "Calling MIGRAD with call limit " << overrideCallLimit << std::endl;
        double plist[1] = {overrideCallLimit};
        int err         = 0;

        if(_useHesseBefore)
            minuit_.mnexcm("HESSE", plist, 1, err);

        minuit_.mnexcm("MIGRAD", plist, 1, err);

        if(_useHesse)
            minuit_.mnexcm("HESSE", plist, 1, err);

        if(_useMinos)
            minuit_.mnexcm("MINOS", plist, 1, err);

        if(_useImprove)
            minuit_.mnexcm("IMPROVE", plist, 1, err);
    } else
        minuit_.Migrad();

    std::cout << GooFit::reset;

    double tmp_value, tmp_error;
    for(Variable &var : minuit_.getVaraibles()) {
        minuit_.GetParameter(var.getFitterIndex(), tmp_value, tmp_error);
        var.setValue(tmp_value);
        var.setError(tmp_error);
    }
}

void FitManagerMinuit1::getMinuitStatus(
    double &fmin, double &fedm, double &errdef, int &npari, int &nparx, int &istat) {
    minuit_.mnstat(fmin, fedm, errdef, npari, nparx, istat);
    std::cout << "mnstat(fmin = " << fmin << ", fedm = " << fedm << ", errdef = " << errdef << ", npari = " << npari
              << ", nparx = " << nparx << ", istat = " << istat << ")" << std::endl;
}
} // namespace GooFit
