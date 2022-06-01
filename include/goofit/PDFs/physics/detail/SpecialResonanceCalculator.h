#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class SpecialResonanceCalculator : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fpcomplex> {
  public:
    // Used to create the cached BW values.
    SpecialResonanceCalculator(int pIdx, unsigned int res_idx);
    void setDalitzIndex(unsigned int id) { dalitz_i = id; }
    void setResonanceIndex(unsigned int id) { resonance_i = id; }
    __device__ auto operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex;

  private:
    unsigned int dalitz_i;
    unsigned int resonance_i;
    unsigned int parameters;
};

} // namespace GooFit
