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

class ThreeGaussResolutionSplice : public MixingTimeResolution {
  public:
    ThreeGaussResolutionSplice(Variable cf,
                               Variable tf,
                               Variable cb,
                               Variable cs,
                               Variable tb,
                               Variable ts,
                               Variable ob,
                               Variable os,
                               std::vector<Variable> knots,
                               std::vector<Variable> a0,
                               std::vector<Variable> a1,
                               std::vector<Variable> a2,
                               std::vector<Variable> a3);
    ~ThreeGaussResolutionSplice() override;

    fptype normalization(
        fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override;

  private:
    std::vector<Variable> m_knots;
    std::vector<Variable> m_a0;
    std::vector<Variable> m_a1;
    std::vector<Variable> m_a2;
    std::vector<Variable> m_a3;
};
} // namespace GooFit
