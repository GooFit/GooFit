#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Cubic spline constructor
class Spline : public ResonancePdf {
  public:
    Spline(std::string name,
           Variable ar,
           Variable ai,
           std::vector<fptype> &HH_bin_limits,
           std::vector<Variable> &pwa_coefs_reals,
           std::vector<Variable> &pwa_coefs_imags,
           unsigned int cyc,
           bool symmDP = false);
    ~Spline() override = default;

    /// Recacluate the CACHE values before running
    __host__ void recalculateCache() const override;
};
} // namespace Resonances

} // namespace GooFit
