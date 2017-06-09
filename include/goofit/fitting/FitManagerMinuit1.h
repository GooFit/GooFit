#pragma once

#include "goofit/GlobalCudaDefines.h"
#include "goofit/Variable.h"
#include <TMinuit.h>

namespace GooFit {

class PdfBase;

class Minuit1 : public TMinuit {
    PdfBase *pdfPointer;
    std::vector<Variable *> vars;

  public:
    Minuit1(PdfBase *pdfPointer);
    /// Fit function for Minuit
    virtual Int_t Eval(Int_t npar,     //< The number of parameters
                       Double_t *grad, //< The derivatives can be stored here if flag is 2 (output)
                       Double_t &fval, //< The value of the function at this point (output)
                       Double_t *par,  //< The input parameters
                       Int_t flag //< This is 1 the first time, 2, for derivatives, and 3 after the fit is finished. It
                                  // is something else if computing.
                       ) override;

    // Get a copy of the list of variables
    std::vector<Variable *> getVaraibles() const { return vars; };
};

class FitManagerMinuit1 {
  public:
    FitManagerMinuit1(PdfBase *dat)
        : minuit_(dat) {}

    void setMaxCalls(double mxc) { overrideCallLimit = mxc; }
    void useHesseBefore(bool use = true) { _useHesseBefore = use; }
    void useHesse(bool use = true) { _useHesse = use; }
    void useMinos(bool use = true) { _useMinos = use; }
    void useImprove(bool use = true) { _useImprove = use; }
    void setVerbosity(int v) { minuit_.SetPrintLevel(v - 1); }

    operator bool() const { return minuit_.GetStatus() == 0; }
    operator int() const { return minuit_.GetStatus(); }

    // This runs the fit
    void fit();

    Minuit1 *getMinuitObject() { return &minuit_; }

    void getMinuitStatus(double &fmin, double &fedm, double &errdef, int &npari, int &nparx, int &istat);

  private:
    double overrideCallLimit{-1};
    bool _useHesseBefore{true};
    bool _useHesse{true};
    bool _useMinos{false};
    bool _useImprove{false};

    Minuit1 minuit_;
};
}
