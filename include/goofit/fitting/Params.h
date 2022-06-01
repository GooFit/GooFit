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

    bool do_record_{false};
    std::vector<std::vector<double>> recorded_;

  public:
    using MnUserParameters::MnUserParameters;

    Params(PdfBase &pdf);

    // Get vector of GooFit Variables.
    std::vector<Variable> GetGooFitParams() const { return vars_; }

    /// Read the values back into GooFit
    void SetGooFitParams(const Minuit2::MnUserParameterState &input);

    /// Get the number of params in the fit
    auto size() const -> size_t { return vars_.size(); };

    /// Make a parameter array with the current variable values
    auto make_minuit_vector() const -> std::vector<double>;

    /// Set from a minuit vector. Optional force_changed to force complete recalculation
    void from_minuit_vector(const std::vector<double> &values, bool force_changed = false);

    /// Set recording of Minuit parameter changes
    void set_record(bool do_record = true) { do_record_ = do_record; }

    /// Get recorded values array
    auto get_recorded() const -> std::vector<std::vector<double>> { return recorded_; }
};
} // namespace GooFit
