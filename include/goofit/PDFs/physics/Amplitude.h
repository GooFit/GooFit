#pragma once

#include <goofit/PDFs/physics/AmpComponent.h>
#include <goofit/Variable.h>
#include <goofit/detail/Complex.h>

#include <vector>

namespace GooFit {

class Lineshape;
class SpinFactor;

class Amplitude : public AmpComponent {
    friend class Amp4Body;
    friend class Amp4Body_TD;

  public:
    Amplitude(std::string uniqueDecayStr,
              Variable ar,
              Variable ai,
              std::vector<Lineshape *> LS,
              std::vector<SpinFactor *> SF,
              unsigned int nPerm = 1);

    auto operator==(const Amplitude &A) const -> bool;

    __host__ void recursiveSetIndices() override;

    auto getSpinFactors() -> std::vector<SpinFactor *> { return _SF; }
    auto getLineShapes() -> std::vector<Lineshape *> { return _LS; }

    // have the lineshape parameters changed since the last evaluation?
    bool lineshapeParametersChanged() const;

  private:
    std::string _uniqueDecayStr;
    Variable _ar;
    Variable _ai;
    std::vector<SpinFactor *> _SF;
    std::vector<Lineshape *> _LS;
    unsigned int _nPerm;
};

} // namespace GooFit
