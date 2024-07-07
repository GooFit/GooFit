#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceCalculator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

SpecialResonanceCalculator::SpecialResonanceCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ auto SpecialResonanceCalculator::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    ParameterContainer pc;

    while(pc.funcIdx < dalitz_i)
        pc.incrementIndex();

    int id_m13 = pc.getObservable(0);
    int id_m23 = pc.getObservable(1);

    fptype m13 = evt[id_m13];
    fptype m23 = evt[id_m23];

    if(!inDalitz2(m13, m23, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return fpcomplex(0.0,0.0);

    if(norm < 1E-10)
        return ret;

    // m12, m23 and m13 stand for the squared invariant masses.
    // Now fixed.
    fptype m12 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - m13 - m23;

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    //printf("s12=%.2f \t s13=%.2f \t s23=%.2f \n",m12, m13, m23);


    ret = getResonanceAmplitude(m13, m23, m12, pc);

  //  printf("Resonance norm = %f \n",sqrt(ret.real()*ret.real() + ret.imag()*ret.imag()));

    return ret/sqrt(norm);
}

} // namespace GooFit
