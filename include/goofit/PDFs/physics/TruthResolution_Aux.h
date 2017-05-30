#pragma once

#include "goofit/PDFs/MixingTimeResolution_Aux.h"

namespace GooFit {


class TruthResolution : public MixingTimeResolution {
public:
    TruthResolution();
    ~TruthResolution();

    virtual fptype normalisation(fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing,
                                 fptype ymixing) const;
    virtual void createParameters(std::vector<unsigned int>& pindices, PdfBase* dis) {
        pindices.push_back(0);
    }
};
} // namespace GooFit

