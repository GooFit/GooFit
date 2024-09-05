#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

namespace GooFit {

namespace Resonances {

/// Nonresonant constructor
class NonRes : public ResonancePdf {
  public:
    NonRes(std::string name, Variable ar, Variable ai);
    ~NonRes() override = default;
};

} // namespace Resonances

} // namespace GooFit
