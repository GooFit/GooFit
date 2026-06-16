#pragma once

#include <vector>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <mcbooster/GTypes.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class AmpCalc_TD : public thrust::unary_function<unsigned int, fpcomplex> {
  public:
    AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx);
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    // void setAmplitudeId(int idx) { _AmpIdx = idx; }
    // void setpIdx(unsigned int pIdx){_parameters = pIdx;}
    // Called via a unary transform over a counting iterator, so it takes the
    // event index directly. (Older Thrust silently let an int convert to a
    // single-element-initialized tuple; CCCL 2.x no longer does.)
    __device__ auto operator()(unsigned int evtNum) const -> fpcomplex;

    __host__ std::vector<unsigned int> getLineshapeIndices(int totalAMP) const;

    __host__ std::vector<unsigned int> getSpinFactorIndices(int totalAMP) const;

  private:
    unsigned int _nPerm;
    unsigned int _AmpIdx;
    unsigned int dalitzFuncId;
};

} // namespace GooFit
