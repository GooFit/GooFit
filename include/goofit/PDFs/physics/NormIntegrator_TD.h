#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class NormIntegrator_TD : public thrust::unary_function<thrust::tuple<int, int, fptype *, fpcomplex *>, fptype> {
  public:
    NormIntegrator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    __device__ thrust::tuple<fptype, fptype, fptype, fptype>
    operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const;

  private:
    unsigned int dalitzFuncId;
};

} // namespace GooFit
