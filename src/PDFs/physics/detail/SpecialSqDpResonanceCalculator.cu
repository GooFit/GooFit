#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/detail/SpecialSqDpResonanceCalculator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/Amp3BodySqDP.h>

namespace GooFit {

SpecialSqDpResonanceCalculator::SpecialSqDpResonanceCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ auto SpecialSqDpResonanceCalculator::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    ParameterContainer pc;

    while(pc.funcIdx < dalitz_i)
        pc.incrementIndex();

    int id_mprime = pc.getObservable(0);
    int id_thetaprime = pc.getObservable(1);

    fptype mprime = RO_CACHE(evt[id_mprime]);
    fptype thetaprime = RO_CACHE(evt[id_thetaprime]);

    if(!inSqDalitz(mprime, thetaprime))
        return ret;

    // mprime, m23 and thetaprime stand for the squared invariant masses.
    // Now fixed.
    fptype m12 = calc_m12(mprime,c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype m13 = calc_m13(thetaprime, m12, c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype s12 = m12*m12;
    fptype s13 = m13*m13;
    fptype s23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - s12 - s13;

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    printf("s12=%.2f \t s13=%.2f \t s23=%.2f \t mp=%.2f \t tp=%.2f \n",s12, s13, s23, mprime, thetaprime);

    ret = getResonanceAmplitude(s12, s13, s23, pc);

    //fptype jacobian = calc_SqDp_Jacobian(mprime, thetaprime, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass);
    //ret*=jacobian;

    return ret;
}

} // namespace GooFit
