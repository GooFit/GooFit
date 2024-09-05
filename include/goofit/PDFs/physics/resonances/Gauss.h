#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

namespace GooFit {

namespace Resonances {

/// Gaussian constructor
class Gauss : public ResonancePdf {
  public:
    Gauss(std::string name,
          Variable ar,
          Variable ai,
          Variable mean,
          Variable sigma,
          unsigned int cyc,
          bool symmDP = false);
    ~Gauss() override = default;
};

} // namespace Resonances

} // namespace GooFit
