#pragma once

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnScan.h>

#include <memory>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/fitting/FCN.h>

namespace GooFit {

class PdfBase;

enum class FitErrors : int { Valid = 0, NotRun = 50, InValid };

class FitManagerMinuit2 {
  public:
    FitManagerMinuit2(PdfBase *dat);

    /// This runs the fit
    auto fit() -> ROOT::Minuit2::FunctionMinimum;

    /// Set the maximum number of calls. 0 for Minuit2 default.
    void setMaxCalls(unsigned int max_calls = 0) { maxfcn_ = max_calls; }

    /// Get a pointer to the params
    auto getParams() -> Params * { return &upar_; }

    /// Get a pointer to the fcn
    auto getFCN() -> FCN * { return &fcn_; }

    ROOT::Minuit2::MnScan getMnScan();

    /// Check to see if fit is valid
    operator bool() const { return retval_ == FitErrors::Valid; }

    /// Return value for program
    operator int() const { return static_cast<int>(retval_); }

    /// Set the fitting verbosity
    void setVerbosity(int value) { verbosity = value; }

    /// Get the fitting verbosity
    auto getVerbosity() const -> int { return verbosity; }

    // Get the minos errors
    std::vector<std::pair<double, double>> getMinosErrors() const { return minos_errors; }

    // Run Minos error calculation
    void setMinos(bool minos_flag = 1) { minos = minos_flag; }

  private:
    Params upar_;
    FCN fcn_;
    unsigned int maxfcn_{0};
    FitErrors retval_{FitErrors::NotRun};
    int verbosity{3};
    std::vector<std::pair<double, double>> minos_errors;
    bool minos{0};
};
} // namespace GooFit
