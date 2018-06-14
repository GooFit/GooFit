#pragma once

#include <thrust/tuple.h>

namespace GooFit {

class LSCalculator_TD : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fpcomplex> {
  public:
    // Used to create the cached BW values.
    LSCalculator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setResonanceId(int idx) { _resonance_i = idx; }
    __device__ fpcomplex operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int _resonance_i{0};
    unsigned int dalitzFuncId;
};

} // namespace GooFit
