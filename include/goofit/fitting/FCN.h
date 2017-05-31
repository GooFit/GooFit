#pragma once

#include <Minuit2/FCNBase.h>

#include <vector>

#include "goofit/fitting/Params.h"

namespace GooFit {

class FCN : public Minuit2::FCNBase {
  protected:
    Params *params_;

  public:
    /// Create an FCN given parameters (PDF reference is inside params)
    FCN(Params &params);

    /// Make a parameter array with the current variable values
    std::vector<double> makePars() const;

    /// Run the fit (used by Minuit2 class)
    double operator()(const std::vector<double> &pars) const override;

    /// produce the FCN value for the current values of the parameters
    double operator()() const;

    /// This value is 0.5 for ll, 1 for chi2
    double Up() const override { return 0.5; }

    /// Get a pointer to the parameters
    Params *GetParams();
};
}
