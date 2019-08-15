#pragma once

#include <goofit/PDFs/physics/AmpNBodyBase.h>
#include <mcbooster/GContainers.h>

namespace GooFit {

class Amp4BodyBase : public AmpNBodyBase {
  public:
    using AmpNBodyBase::AmpNBodyBase;

    /// Make the accept/reject flags from weights
    __host__ void
    fillMCFlags(mcbooster::BoolVector_d &flags, const mcbooster::RealVector_d &weights, unsigned int numEvents);
};

} // namespace GooFit
