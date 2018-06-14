#pragma once

#include <thrust/tuple.h>

namespace GooFit {

class AmpCalc_TD : public thrust::unary_function<unsigned int, fpcomplex> {
  public:
    AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx);
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    // void setAmplitudeId(int idx) { _AmpIdx = idx; }
    // void setpIdx(unsigned int pIdx){_parameters = pIdx;}
    __device__ fpcomplex operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int _nPerm;
    unsigned int _AmpIdx;
    unsigned int dalitzFuncId;
};

} // namespace GooFit
