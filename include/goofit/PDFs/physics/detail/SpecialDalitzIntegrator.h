#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

class SpecialDalitzIntegrator : public thrust::unary_function<thrust::tuple<int, fptype *, int>, ThreeComplex> {
  public:
    SpecialDalitzIntegrator(int pIdx, unsigned int ri, unsigned int rj);
    void setTddpIndex(unsigned int id) { tddp = id; }
    void setResonanceIndex(unsigned int id) { resonance_i = id; }
    void setEfficiencyIndex(unsigned int id) { resonance_j = id; }
    __device__ auto operator()(thrust::tuple<int, fptype *, int> t) const -> ThreeComplex;

  private:
    unsigned int tddp;
    unsigned int resonance_i;
    unsigned int resonance_j;
    unsigned int parameters;
};

} // namespace GooFit
