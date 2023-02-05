#include <goofit/PDFs/physics/detail/SpecialIncoherentResonanceCalculator.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

SpecialIncoherentResonanceCalculator::SpecialIncoherentResonanceCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ auto SpecialIncoherentResonanceCalculator::operator()(thrust::tuple<int, fptype *, int> t) const
    -> fpcomplex {
    // Returns the BW, or other resonance function, for a specific resonance.
    // Is special because the value is expected to change slowly, so it's
    // useful to cache the result.
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    // unsigned int *indices = paramIndices + parameters; // Jump to TDDP position within parameters array
    ParameterContainer pc;

    while(pc.funcIdx < incoherentSum)
        pc.incrementIndex();

    int id_m12 = pc.getObservable(0);
    int id_m13 = pc.getObservable(1);

    fptype m12 = RO_CACHE(evt[id_m12]);
    fptype m13 = RO_CACHE(evt[id_m13]);

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return {0.0, 0.0};

    fptype m23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - m12 - m13;

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    // int parameter_i
    //    = parIndexFromResIndex_incoherent(resonance_i); // Find position of this resonance relative to TDDP start
    // unsigned int functn_i       = indices[parameter_i + 2];
    // unsigned int params_i       = indices[parameter_i + 3];
    fpcomplex ret = getResonanceAmplitude(m12, m13, m23, pc);

    return ret;
}

} // namespace GooFit
