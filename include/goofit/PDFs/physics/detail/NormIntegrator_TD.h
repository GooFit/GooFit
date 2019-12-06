#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>
#include <mcbooster/GContainers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
namespace GooFit {

class NormIntegrator_TD :
  // public thrust::unary_function<thrust::tuple<int, int, fptype *, fpcomplex *>, fptype> 
  public thrust::unary_function< thrust::tuple<int, int, fptype *, thrust::complex<fptype>*,
  mcbooster::GReal_t, mcbooster::GReal_t>,
  thrust::tuple<fptype, fptype, fptype, fptype>>{
  public:
    NormIntegrator_TD(bool SpecInt);
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    //__device__ thrust::tuple<fptype, fptype, fptype, fptype>
    //operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const;
    __device__ thrust::tuple<fptype, fptype, fptype, fptype>
      operator()(thrust::tuple<int, int, fptype *, thrust::complex<fptype>*,
		 mcbooster::GReal_t, mcbooster::GReal_t> t) const;
  private:
    unsigned int dalitzFuncId;
    bool _SpecInt;
};

} // namespace GooFit
