#pragma once

#include <goofit/Variable.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>

#include <vector>

namespace Minuit2 = ROOT::Minuit2;

namespace GooFit {

class PdfBase;
class FCN;

class Params : public Minuit2::MnUserParameters {
    friend FCN;

  protected:
    std::vector<Variable> vars_;
    PdfBase *pdf_;
    size_t num_;

  public:
    using MnUserParameters::MnUserParameters;

    Params(PdfBase &pdf);

    /// Read the values back into GooFit
    void SetGooFitParams(const Minuit2::MnUserParameterState &input);

    /// Get the number of params in the fit
    size_t size() const { return vars_.size(); };

    /// Make a parameter array with the current variable values
    std::vector<double> make_minuit_vector() const;

    /// Set from a minuit vector. Optional force_changed to force complete recalculation
    void from_minuit_vector(const std::vector<double> &values, bool force_changed = false);
};
} // namespace GooFit
