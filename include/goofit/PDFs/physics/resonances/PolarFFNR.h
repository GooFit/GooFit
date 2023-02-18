#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Relativistic Breit-Wigner
class PolarFFNR : public ResonancePdf {
  public:
    PolarFFNR(std::string name,
        Variable ar,
        Variable ai,
        Variable lambda,
        unsigned int cyc,
        bool sym = false);
    ~PolarFFNR() override = default;
};

} // namespace Resonances

} // namespace GooFit
