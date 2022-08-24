#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Rho-omega mixing
class RhoOmegaMix : public ResonancePdf {
  public:
    RhoOmegaMix(std::string name,
                Variable ar,
                Variable ai,
                Variable omega_mass,
                Variable omega_width,
                Variable rho_mass,
                Variable rho_width,
                Variable real,
                Variable imag,
                Variable delta,
                unsigned int sp,
                unsigned int cyc,
                bool sym,
                bool bachPframe=false,
                bool ignoreMom=false,
                bool ignoreBW=false);

    ~RhoOmegaMix() override = default;
};

} // namespace Resonances

} // namespace GooFit
