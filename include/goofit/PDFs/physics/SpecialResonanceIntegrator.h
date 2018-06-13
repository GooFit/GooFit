#pragma once

#include <thrust/tuple.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

class SpecialResonanceIntegrator : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fpcomplex> {
  public:
    // Class used to calculate integrals of terms BW_i * BW_j^*.
    SpecialResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj);
    void setDalitzIndex(unsigned int id) { dalitz_i = id; }
    void setResonanceIndex(unsigned int id) { resonance_i = id; }
    void setEfficiencyIndex(unsigned int id) { resonance_j = id; }
    __device__ fpcomplex operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int dalitz_i;
    unsigned int resonance_i;
    unsigned int resonance_j;
    unsigned int parameters;
};

} // namepace GooFit
