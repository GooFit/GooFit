#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <mcbooster/GTypes.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class NormSpinCalculator_TD
    : public thrust::unary_function<thrust::tuple<fptype, fptype, fptype, fptype, fptype>, fptype> {
  public:
    // Used to create the cached BW values.
    NormSpinCalculator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setSpinFactorId(int idx) { _spinfactor_i = idx; }
    __device__ auto operator()(thrust::tuple<fptype, fptype, fptype, fptype, fptype> t) const -> fptype;

  private:
    unsigned int _spinfactor_i{0};
    unsigned int dalitzFuncId;
};

} // namespace GooFit
