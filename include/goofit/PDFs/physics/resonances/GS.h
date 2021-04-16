#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Gounaris-Sakurai
class GS : public ResonancePdf {
  public:
    GS(std::string name,
       Variable ar,
       Variable ai,
       Variable mass,
       Variable width,
       unsigned int sp,
       unsigned int cyc,
       bool norm = true,
       bool sym  = false);
    ~GS() override = default;
};

} // namespace Resonances

} // namespace GooFit
