#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Relativistic Breit-Wigner
class RBW : public ResonancePdf {
  public:
    RBW(std::string name,
        Variable ar,
        Variable ai,
        Variable mass,
        Variable width,
        unsigned int sp,
        unsigned int cyc,
        bool norm = true,
        bool sym  = false);
    ~RBW() override = default;
};

} // namespace Resonances

} // namespace GooFit
