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
    __device__ auto operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex;

    __host__ std::vector<unsigned int> getLineshapeIndices(int totalAMP) const;

    __host__ std::vector<unsigned int> getSpinFactorIndices(int totalAMP) const;

  private:
    unsigned int _nPerm;
    unsigned int _AmpIdx;
    unsigned int dalitzFuncId;
};

} // namespace GooFit
