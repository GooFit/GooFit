#pragma once

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnPlot.h>

#include <memory>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/fitting/FCN.h>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using namespace Eigen;

namespace GooFit {

class PdfBase;

enum class FitErrors : int { Valid = 0, NotRun = 50, InValid };

class FitManagerMinuit2 {
  public:
    FitManagerMinuit2(PdfBase *dat);

    /// This runs the fit
    auto fit() -> ROOT::Minuit2::FunctionMinimum;

    // This runs the scan
    auto scan(unsigned int par_i, unsigned int maxsteps = 41, fptype low = 0., fptype high = 0.)
        -> std::vector<std::pair<fptype, fptype>>;

    auto getPdf() const -> PdfBase * { return pdfPointer; }

    /// Set the maximum number of calls. 0 for Minuit2 default.
    void setMaxCalls(unsigned int max_calls = 0) { maxfcn_ = max_calls; }

    void setStrategy(unsigned int str) { strategy = str; }

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

    // Set minuit tolerance
    void setTolerance(fptype tolerance) { tolerance_ = tolerance; }

    auto getTolerance() const -> fptype { return tolerance_; }

    // Convert real and imag coefs to mag and phase.
    void printCovMat();
    auto dpda(fptype, fptype) -> fptype;
    auto dpdb(fptype, fptype) -> fptype;
    auto dmda(fptype, fptype) -> fptype;
    auto dmdb(fptype, fptype) -> fptype;
    std::vector<std::vector<fptype>> printParams();
    void printOriginalParams();
    void setRandMinuitValues(size_t nSamples);
    void loadSample(size_t iSample);

    // Get the minos errors
    std::vector<std::pair<double, double>> getMinosErrors() const { return minos_errors; }

    // Run Minos error calculation
    void setMinos(bool minos_flag = 1) { minos = minos_flag; }

  private:
    unsigned int strategy = 1;
    Params upar_;
    PdfBase *pdfPointer;
    FCN fcn_;
    unsigned int maxfcn_{0};
    FitErrors retval_{FitErrors::NotRun};
    int verbosity{3};

    fptype tolerance_{0.1};
    Minuit2::MnUserCovariance matCov;
    MatrixXd *sqrtCov;
    std::vector<VectorXd> samples;

    std::vector<std::pair<double, double>> minos_errors;
    bool minos{0};
};
} // namespace GooFit
