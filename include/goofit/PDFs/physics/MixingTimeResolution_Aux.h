#pragma once

#include "goofit/GlobalCudaDefines.h"
#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

typedef fptype (*device_resfunction_ptr)(
    fptype, fptype, fptype, fptype, fptype, fptype, fptype, fptype, fptype, fptype *, unsigned int *);
typedef fptype (*device_calc_tau_fcn_ptr)(fptype, fptype, fptype, fptype, fptype, fptype, fptype);

class MixingTimeResolution {
  public:
    MixingTimeResolution();
    ~MixingTimeResolution();

    void initIndex(void *dev_fcn_ptr = host_fcn_ptr);

    virtual fptype
    normalisation(fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const = 0;
    virtual void createParameters(std::vector<unsigned int> &pindices, PdfBase *dis) = 0;
    int getDeviceFunction() const { return resFunctionIdx; }
    int getCalcTauIdx() const { return resCalcTauFcnIdx; }
    void setCalcTauIdx(int idx) { resCalcTauFcnIdx = idx; }

  private:
    int resFunctionIdx;
    int resCalcTauFcnIdx;
};
} // namespace GooFit
