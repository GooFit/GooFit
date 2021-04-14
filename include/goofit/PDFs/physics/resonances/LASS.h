#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// LASS
class LASS : public ResonancePdf {
  public:
    LASS(std::string name,
         Variable ar,
         Variable ai,
         Variable mass,
         Variable width,
         Variable _a,
         Variable _r,
         Variable _R,
         Variable _phiR,
         Variable _B,
         Variable _phiB,
         unsigned int sp,
         unsigned int cyc,
         bool norm = true);
    ~LASS() override = default;
};

} // namespace Resonances

} // namespace GooFit
