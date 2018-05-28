#pragma once

#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>

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

    fptype normalisation(
        fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override;
    void createParameters(PdfBase *dis) override {}
};
} // namespace GooFit
