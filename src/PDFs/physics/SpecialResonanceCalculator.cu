#include <goofit/PDFs/physics/SpecialResonanceCalculator.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

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

        fptype m12 = evt[id_m12];
        fptype m13 = evt[id_m13];

        fptype motherMass = c_motherMass; // pc.constants[pc.constantIdx + 4];
        fptype daug1Mass  = c_daug1Mass;  // pc.constants[pc.constantIdx + 5];
        fptype daug2Mass  = c_daug2Mass;  // pc.constants[pc.constantIdx + 6];
        fptype daug3Mass  = c_daug3Mass;  // pc.constants[pc.constantIdx + 7];

        if(!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass))
            return ret;

        fptype m23
            = motherMass * motherMass + daug1Mass * daug1Mass + daug2Mass * daug2Mass + daug3Mass * daug3Mass - m12 - m13;

        while(pc.funcIdx < resonance_i)
            pc.incrementIndex();

        ret = getResonanceAmplitude(m12, m13, m23, pc);

        return ret;
    }

}
