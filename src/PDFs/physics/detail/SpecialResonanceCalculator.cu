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
    fpcomplex ret_mt;

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

    // m12, m23 and m13 stand for the squared invariant masses.
    // Now fixed.
    fptype m23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                      + c_daug3Mass * c_daug3Mass - m12 - m13;

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    fptype mistag = 0.005;

    ParameterContainer tmp = pc;
    ret                    = getResonanceAmplitude(m12, m13, m23, tmp);
    tmp                    = pc;
    ret_mt                 = getResonanceAmplitude(m13, m12, m23, tmp);

    ret *= (1-mistag);
    ret_mt *= mistag;
    ret += ret_mt;

    return ret;
}

} // namespace GooFit
