#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3Body_TD.h>
#include <goofit/PDFs/physics/detail/SpecialWaveCalculator.h>

#include <thrust/transform_reduce.h>

namespace GooFit {

SpecialWaveCalculator::SpecialWaveCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ auto SpecialWaveCalculator::operator()(thrust::tuple<int, fptype *, int> t) const -> WaveHolder_s {
    // Calculates the BW values for a specific resonance.
    // The 'A' wave stores the value at each point, the 'B'
    // at the opposite (reversed) point.

    WaveHolder_s ret;
    ret.ai_real = 0.0;
    ret.ai_imag = 0.0;
    ret.bi_real = 0.0;
    ret.bi_imag = 0.0;

    int evtNum  = thrust::get<0>(t);
    int evtSize = thrust::get<2>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * evtSize);

    ParameterContainer pc;

    // increment until we are at tddp index
    while(pc.funcIdx < tddp)
        pc.incrementIndex();

    int id_m12 = pc.getObservable(2);
    int id_m13 = pc.getObservable(3);

    // Read these values as tddp.
    fptype m12 = RO_CACHE(evt[id_m12]);
    fptype m13 = RO_CACHE(evt[id_m13]);

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return ret;

    fptype m23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - m12 - m13;

    // int parameter_i       = parIndexFromResIndex(resonance_i); // Find position of this resonance relative to TDDP
    // start  unsigned int functn_i = indices[parameter_i + 2];  unsigned int params_i = indices[parameter_i + 3];

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    ParameterContainer tmp = pc;
    fpcomplex ai           = getResonanceAmplitude(m12, m13, m23, tmp);
    tmp                    = pc;
    fpcomplex bi           = getResonanceAmplitude(m13, m12, m23, tmp);

    // printf("Amplitudes %f, %f => (%f %f) (%f %f)\n", m12, m13, ai.real, ai.imag, bi.real, bi.imag);

    ret.ai_real = ai.real();
    ret.ai_imag = ai.imag();
    ret.bi_real = bi.real();
    ret.bi_imag = bi.imag();

    return ret;
}

} // namespace GooFit
