#pragma once

#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>

namespace GooFit {

/**
A resolution function consisting of a
sum of three Gaussians, referred to as the 'core', 'tail', and
'outlier' components. The constructor takes the core and tail
fractions (the outlier fraction is 1 minus the other two), core mean
and width, tail mean and width, and outlier mean and width. Notice
that this is a resolution function, so the full probability is found
by convolving Gaussians with the equation from GooFit::Amp3BodyT, and this runs
to a page or so of algebra involving error functions. It is beyond
the scope of this documentation.
**/

class ThreeGaussResolution : public MixingTimeResolution {
  public:
    ThreeGaussResolution(
        Variable cf, Variable tf, Variable cb, Variable cs, Variable tb, Variable ts, Variable ob, Variable os);
    ~ThreeGaussResolution() override;

    fptype normalization(
        fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override;
    void createParameters(PdfBase *dis) override;

  private:
    Variable coreFraction;
    Variable tailFraction;
    Variable coreBias;
    Variable coreScaleFactor;
    Variable tailBias;
    Variable tailScaleFactor;
    Variable outBias;
    Variable outScaleFactor;
};
} // namespace GooFit
