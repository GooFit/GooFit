#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

#include <mcbooster/Generate.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/EvaluateArray.h>


namespace GooFit {

class NormIntegrator_TD : public thrust::unary_function<thrust::tuple<int, int, fptype *, fpcomplex *>, fptype> {
  public:
    NormIntegrator_TD(unsigned int CacheIdx);
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    __device__ auto operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const
        -> thrust::tuple<fptype, fptype, fptype, fptype>;

  private:
    unsigned int dalitzFuncId;
    unsigned int _CacheIdx;
};

} // namespace GooFit
