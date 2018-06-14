#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Gaussian constructor
class Gauss : public ResonancePdf {
  public:
    Gauss(std::string name, Variable ar, Variable ai, Variable mean, Variable sigma, unsigned int cyc);
    ~Gauss() override = default;
};

} // namespace Resonances

} // namespace GooFit
