#pragma once

#include <thrust/tuple.h>

namespace GooFit {

class NormIntegrator : public thrust::unary_function<thrust::tuple<int, int, fptype *, fpcomplex *>, fptype> {
  public:
    NormIntegrator();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    __device__ fptype operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const;

  private:
    unsigned int dalitzFuncId;
};

} // namespace GooFit
