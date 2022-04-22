#pragma once

#include <Minuit2/FunctionMinimum.h>

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

    auto getPdf() const -> PdfBase* { return pdfPointer;}

    /// Set the maximum number of calls. 0 for Minuit2 default.
    void setMaxCalls(unsigned int max_calls = 0) { maxfcn_ = max_calls; }

    /// Get a pointer to the params
    auto getParams() -> Params * { return &upar_; }

    /// Get a pointer to the fcn
    auto getFCN() -> FCN * { return &fcn_; }

    /// Check to see if fit is valid
    operator bool() const { return retval_ == FitErrors::Valid; }

    /// Return value for program
    operator int() const { return static_cast<int>(retval_); }

    /// Set the fitting verbosity
    void setVerbosity(int value) { verbosity = value; }

    /// Get the fitting verbosity
    auto getVerbosity() const -> int { return verbosity; }

    //Set minuit tolerance
    void setTolerance(fptype tolerance){ tolerance_=tolerance;}

    auto getTolerance() const -> fptype{return tolerance_;}

    //Convert real and imag coefs to mag and phase.
    void printCovMat();
    auto dpda(fptype, fptype) -> fptype;
    auto dpdb(fptype, fptype) -> fptype;
    auto dmda(fptype, fptype) -> fptype;
    auto dmdb(fptype, fptype) -> fptype;
    std::vector <std::vector<fptype>> printParams();
    void printOriginalParams();
    void setRandMinuitValues (size_t nSamples);
    void loadSample(size_t iSample);

  private:
    Params upar_;
    PdfBase *pdfPointer;
    FCN fcn_;
    unsigned int maxfcn_{0};
    FitErrors retval_{FitErrors::NotRun};
    int verbosity{3};
    fptype tolerance_{0.1};
    Minuit2::MnUserCovariance matCov;
    MatrixXd* sqrtCov;
    std::vector<VectorXd> samples;
};
} // namespace GooFit
