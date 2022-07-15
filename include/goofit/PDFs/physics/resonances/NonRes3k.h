#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Nonresonant constructor
class NonRes3k : public ResonancePdf {
  public:
    NonRes3k(std::string name, Variable ar, Variable ai, Variable alpha);
    ~NonRes3k() override = default;
};

} // namespace Resonances

} // namespace GooFit
