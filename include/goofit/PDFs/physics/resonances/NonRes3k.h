#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

namespace GooFit {

namespace Resonances {

/// Nonresonant constructor
class NonRes3k : public ResonancePdf {
  public:
    NonRes3k(std::string name, Variable ar, Variable ai, Variable alpha, Variable beta);
    ~NonRes3k() override = default;
};

} // namespace Resonances

} // namespace GooFit
