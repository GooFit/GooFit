#include <goofit/PDFs/physics/detail/LSCalculator_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

LSCalculator_TD::LSCalculator_TD() = default;

__device__ fpcomplex LSCalculator_TD::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;

    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    // unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    // int parameter_i
    //    = 12 + (2 * indices[6]) + (_resonance_i * 2); // Find position of this resonance relative to DALITZPLOT start
    // unsigned int functn_i = indices[parameter_i];
    // unsigned int params_i = indices[parameter_i + 1];
    // unsigned int pair     = (paramIndices + params_i)[5];

    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    int id_m12   = pc.getObservable(0);
    int id_m34   = pc.getObservable(1);
    int id_cos12 = pc.getObservable(2);
    int id_cos34 = pc.getObservable(3);
    int id_phi   = pc.getObservable(4);

    fptype m12   = evt[id_m12];
    fptype m34   = evt[id_m34];
    fptype cos12 = evt[id_cos12];
    fptype cos34 = evt[id_cos34];
    fptype phi   = evt[id_phi];

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);

    while(pc.funcIdx < _resonance_i)
        pc.incrementIndex();

    unsigned int pair = pc.getConstant(0);

    if(pair < 2) {
        fptype mres = pair == 0 ? m12 : m34;
        fptype d1   = pair == 0 ? m1 : m3;
        fptype d2   = pair == 0 ? m2 : m4;
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
        // printf("LS_nt %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    } else {
        fptype vecs[16];
        get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
        // printf("LS %i: mass:%f, %f i%f\n",_resonance_i, mres, ret.real, ret.imag );
    }

    // printf("LS %i: %.7g, %.7g, %.7g, %.7g, %.7g \n",evtNum, m12, m34, cos12, cos34, phi );

    // if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
    // printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    // printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    // printf("BW_%i : %f %f\n", _resonance_i, ret.real, ret.imag);

    return ret;
}

} // namespace GooFit
