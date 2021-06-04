#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <mcbooster/GTypes.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class NormLSCalculator_TD
    : public thrust::unary_function<
          thrust::
              tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t>,
          fpcomplex> {
  public:
    // Used to create the cached BW values.
    NormLSCalculator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setResonanceId(int idx) { _resonance_i = idx; }
    __device__ auto operator()(
        thrust::
            tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
        const -> fpcomplex;

  private:
    unsigned int _resonance_i{0};
    unsigned int dalitzFuncId;
};

} // namespace GooFit
