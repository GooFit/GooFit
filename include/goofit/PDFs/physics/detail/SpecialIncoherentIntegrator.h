#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <thrust/functional.h>

namespace GooFit {

class SpecialIncoherentIntegrator : public thrust::unary_function<thrust::tuple<int, fptype *>, fptype> {
  public:
    SpecialIncoherentIntegrator(int pIdx, unsigned int ri);
    void setIncoherentIndex(const unsigned int idx) { incoherentSum = idx; }
    void setEfficiencyIndex(const unsigned int eff) { efficiency = eff; }
    void setResonanceIndex(const unsigned int res) { resonance_i = res; }
    __device__ auto operator()(thrust::tuple<int, fptype *> t) const -> fptype;

  private:
    unsigned int incoherentSum;
    unsigned int efficiency;
    unsigned int resonance_i;
    unsigned int parameters;
};

} // namespace GooFit
