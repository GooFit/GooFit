#pragma once

#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

/**
A resolution function consisting of a
sum of three Gaussians, referred to as the 'core', 'tail', and
'outlier' components. The constructor takes the core and tail
fractions (the outlier fraction is 1 minus the other two), core mean
and width, tail mean and width, and outlier mean and width. Notice
that this is a resolution function, so the full probability is found
by convolving Gaussians with the equation from GooFit::Amp3Body_TD, and this runs
to a page or so of algebra involving error functions. It is beyond
the scope of this documentation.
**/

class ThreeGaussResolutionExt : public MixingTimeResolution {
  public:
    ThreeGaussResolutionExt(Variable cf,
                            Variable tf,
                            Variable cb,
                            Variable cs,
                            Variable tb,
                            Variable ts,
                            Variable ob,
                            Variable os,
                            Variable sb_low,
                            Variable sb_high,
                            Variable Tthres,
                            Variable constantC);
    ~ThreeGaussResolutionExt() override;

    fptype normalization(
        fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override;

  private:
    Variable selectionBias_low;
    Variable selectionBias_high;
    Variable mTthreshold;
    Variable mConstantC;
};
} // namespace GooFit
