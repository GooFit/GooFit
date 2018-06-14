#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// FLATTE constructor
class FLATTE : public ResonancePdf {
  public:
    FLATTE(std::string name,
           Variable ar,
           Variable ai,
           Variable mean,
           Variable g1,
           Variable rg2og1,
           unsigned int cyc,
           bool symmDP);
    ~FLATTE() override = default;
};

} // namespace Resonances

} // namespace GooFit
