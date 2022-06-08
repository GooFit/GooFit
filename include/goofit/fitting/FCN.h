#pragma once

#include <Minuit2/FCNBase.h>

#include <vector>

#include <goofit/fitting/Params.h>

namespace GooFit {

class FCN : public Minuit2::FCNBase {
  protected:
    Params *params_;

  public:
    /// Create an FCN given parameters (PDF reference is inside params)
    FCN(Params &params);

    /// Run the fit (used by Minuit2 class)
    auto operator()(const std::vector<double> &pars) const -> double override;

    /// produce the FCN value for the current values of the parameters
    auto operator()() const -> double;

    /// This value is 0.5 for ll, 1 for chi2
    auto Up() const -> double override { return 0.5; }

    /// Get a pointer to the parameters
    auto GetParams() -> Params *;
};
} // namespace GooFit
