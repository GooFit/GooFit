#pragma once

#include "goofit/PDFs/physics/MixingTimeResolution_Aux.h"

namespace GooFit {

class ThreeGaussResolution : public MixingTimeResolution {
  public:
    ThreeGaussResolution(
        Variable *cf, Variable *tf, Variable *cb, Variable *cs, Variable *tb, Variable *ts, Variable *ob, Variable *os);
    ~ThreeGaussResolution();

    fptype normalisation(
        fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override;
    void createParameters(std::vector<unsigned int> &pindices, PdfBase *dis) override;

  private:
    Variable *coreFraction;
    Variable *tailFraction;
    Variable *coreBias;
    Variable *coreScaleFactor;
    Variable *tailBias;
    Variable *tailScaleFactor;
    Variable *outBias;
    Variable *outScaleFactor;
};
} // namespace GooFit
