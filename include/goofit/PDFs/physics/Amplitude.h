#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/Variable.h>
#include <goofit/detail/Complex.h>

#include <vector>

namespace GooFit {

class Lineshape;
class SpinFactor;

class Amplitude : public GooPdf {
    friend class Amp4Body;
    friend class Amp4BodyT;

  public:
    Amplitude(std::string uniqueDecayStr,
              Variable ar,
              Variable ai,
              std::vector<Lineshape *> LS,
              std::vector<SpinFactor *> SF,
              unsigned int nPerm = 1);

    bool operator==(const Amplitude &A) const;

    __host__ void recursiveSetIndices() override;

    std::vector<SpinFactor *> getSpinFactors() { return _SF; }
    std::vector<Lineshape *> getLineShapes() { return _LS; }

  private:
    std::string _uniqueDecayStr;
    Variable _ar;
    Variable _ai;
    std::vector<SpinFactor *> _SF;
    std::vector<Lineshape *> _LS;
    unsigned int _nPerm;
};

} // namespace GooFit
