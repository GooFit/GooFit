#pragma once

#include <vector>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/utilities/Cache_SF_LS_TD_EntryFinder.h>

#include <mcbooster/GTypes.h>

#include <thrust/functional.h>
#include <thrust/tuple.h>

namespace GooFit {

class AmpCalc_TD : public thrust::unary_function<unsigned int, fpcomplex> {
  public:
    AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx);

    void setDalitzId(int idx) { _entryFinder.setDalitzId(idx); }

    __device__ auto operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex;

    __host__ std::vector<unsigned int> getLineshapeIndices(int totalAMP) const;

    __host__ std::vector<unsigned int> getSpinFactorIndices(int totalAMP) const;

    thrust::host_vector<fpcomplex>
    debugLS(unsigned int lsNum, const thrust::device_vector<unsigned int> &evtNums) const;

    thrust::host_vector<fptype> debugSF(unsigned int sfNum, const thrust::device_vector<unsigned int> &evtNums) const;

  private:
    Cache_SF_LS_TD_EntryFinder _entryFinder;
};

} // namespace GooFit
