#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

typedef fptype (*device_resfunction_ptr)(
    fptype, fptype, fptype, fptype, fptype, fptype, fptype, fptype, fptype, ParameterContainer &pc);
typedef fptype (*device_calc_tau_fcn_ptr)(fptype, fptype, fptype, fptype, fptype, fptype, fptype);

/**
The abstract base class of
GooFit::TruthResolution and GooFit::ThreeGaussResolution.
Represents a parametrization of the time resolution.
**/

class MixingTimeResolution : public GooPdf {
  public:
    template <class... Args>
    MixingTimeResolution(std::string pdf_name, Args &&...args)
        : GooPdf(pdf_name, "mixing_resolution", std::forward<Args>(args)...) {}

    ~MixingTimeResolution() override;

    void initIndex(void *dev_fcn_ptr = host_fcn_ptr);

    virtual auto
    normalization(fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const
        -> fptype = 0;

    auto getDeviceFunction() const -> int { return resFunctionIdx; }
    auto getCalcTauIdx() const -> int { return resCalcTauFcnIdx; }
    void setCalcTauIdx(int idx) { resCalcTauFcnIdx = idx; }

  private:
    int resFunctionIdx;
    int resCalcTauFcnIdx;
};
} // namespace GooFit
