#pragma once

#include "goofit/PDFs/physics/MixingTimeResolution_Aux.h"

namespace GooFit {

class TruthResolution : public MixingTimeResolution {
  public:
    TruthResolution();
    ~TruthResolution();

    fptype normalisation(
        fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override;
    void createParameters(std::vector<unsigned int> &pindices, PdfBase *dis) override { pindices.push_back(0); }
};
} // namespace GooFit
