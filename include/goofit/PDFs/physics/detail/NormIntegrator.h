#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class NormIntegrator : public thrust::unary_function<thrust::tuple<int, int, fptype *, fpcomplex *>, fptype> {
  public:
    NormIntegrator();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    __device__ auto operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const -> fptype;

  private:
    unsigned int dalitzFuncId;
};

} // namespace GooFit
