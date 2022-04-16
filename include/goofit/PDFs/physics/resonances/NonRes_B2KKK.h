#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// NonRes_B2KKKonant constructor
class NonRes_B2KKK : public ResonancePdf {
  public:
    NonRes_B2KKK(std::string name, Variable ar, Variable ai, Variable alpha, Variable beta);
    ~NonRes_B2KKK() override = default;
};

} // namespace Resonances

} // namespace GooFit
