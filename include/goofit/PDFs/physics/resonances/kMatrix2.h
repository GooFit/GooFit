#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Relativistic Breit-Wigner
class kMatrix2 : public ResonancePdf {
  public:
    kMatrix2(std::string name,
        Variable ar,
        Variable ai,
        Variable mass,
        Variable width,
        unsigned int sp,
        unsigned int cyc,
        bool norm = true,
        bool sym  = false);
    ~kMatrix2() override = default;
};

} // namespace Resonances

} // namespace GooFit
