#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>

namespace GooFit {

class SpecialIncoherentResonanceCalculator
    : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fpcomplex> {
  public:
    SpecialIncoherentResonanceCalculator(int pIdx, unsigned int res_idx);
    void setIncoherentIndex(const unsigned int idx) { incoherentSum = idx; }
    void setResonanceIndex(const unsigned int res) { resonance_i = res; }
    __device__ auto operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex;

  private:
    unsigned int incoherentSum;
    unsigned int resonance_i;
    unsigned int parameters;
};

} // namespace GooFit
