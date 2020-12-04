#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceCalculator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

SpecialResonanceCalculator::SpecialResonanceCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ fpcomplex SpecialResonanceCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    ParameterContainer pc;

    while(pc.funcIdx < dalitz_i)
        pc.incrementIndex();

    int id_m12 = pc.getObservable(0);
    int id_m13 = pc.getObservable(1);

    fptype m12 = RO_CACHE(evt[id_m12]);
    fptype m13 = RO_CACHE(evt[id_m13]);

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return ret;

    // This was being calculated with invariant masses instead of squared invariante masses and without taking the sqrt.
    // Now fixed.
    fptype m23 = sqrt(c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                      + c_daug3Mass * c_daug3Mass - m12 * m12 - m13 * m13);

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    ret = getResonanceAmplitude(m12, m13, m23, pc);

    return ret;
}

} // namespace GooFit
