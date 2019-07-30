#pragma once

#include <goofit/PDFs/physics/AmpNBodyBase.h>
#include <mcbooster/GContainers.h>

namespace GooFit {

class Amp4BodyBase : public AmpNBodyBase {
  public:
    using AmpNBodyBase::AmpNBodyBase;

    /// Make the accept/reject flags from weights
    __host__ mcbooster::BoolVector_d makeMCFlags(const mcbooster::RealVector_d &weights, unsigned int numEvents);
};

} // namespace GooFit
