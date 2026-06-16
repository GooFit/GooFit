#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <mcbooster/GTypes.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class AmpCalc : public thrust::unary_function<unsigned int, fpcomplex> {
  public:
    AmpCalc(unsigned int nPerm, unsigned int ampIdx);
    // void setpIdx(unsigned int pIdx){_parameters = pIdx;}
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setAmplitudeId(int idx) { _AmpIdx = idx; }
    // Called via a unary transform over a counting iterator, so it takes the
    // event index directly. (Older Thrust silently let an int convert to a
    // single-element-initialized tuple; CCCL 2.x no longer does.)
    __device__ auto operator()(unsigned int evtNum) const -> fpcomplex;

  private:
    unsigned int dalitzFuncId;
    unsigned int _nPerm;
    unsigned int _AmpIdx;
};

} // namespace GooFit
