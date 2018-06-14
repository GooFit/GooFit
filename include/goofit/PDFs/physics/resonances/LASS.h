#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// LASS
class LASS : public ResonancePdf {
  public:
    LASS(std::string name, Variable ar, Variable ai, Variable mass, Variable width, unsigned int sp, unsigned int cyc);
    ~LASS() override = default;
};

} // namespace Resonances

} // namespace GooFit
