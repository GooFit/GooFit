#pragma once

#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

/**
The simplest possible resolution function, a
simple delta spike at zero - i.e., time is always measured
perfectly. The constructor takes no arguments at all!
**/

class TruthResolution : public MixingTimeResolution {
  public:
    TruthResolution();
    ~TruthResolution() override;

    auto normalization(fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const
        -> fptype override;
};
} // namespace GooFit
