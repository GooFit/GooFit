#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class NormIntegrator_TD_Weighted
    : public thrust::unary_function<thrust::tuple<int, int, fptype *, fptype *, fptype *, fptype *, fpcomplex *>,
                                    fptype> {
  public:
    NormIntegrator_TD_Weighted();

    void setDalitzId(int idx) { dalitzFuncId = idx; }

    __device__ fptype operator()(thrust::tuple<int, int, fptype *, fptype *, fptype *, fptype *, fpcomplex *> t) const;

  private:
    unsigned int dalitzFuncId;
};

} // namespace GooFit
