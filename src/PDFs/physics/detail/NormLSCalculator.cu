#include <goofit/PDFs/physics/detail/NormLSCalculator.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

NormLSCalculator::NormLSCalculator() = default;

__device__ auto NormLSCalculator::operator()(
    thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const -> fpcomplex {
    // Calculates the BW values for a specific resonance.
    fpcomplex ret;

    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    fptype m12   = (thrust::get<0>(t));
    fptype m34   = (thrust::get<1>(t));
    fptype cos12 = (thrust::get<2>(t));
    fptype cos34 = (thrust::get<3>(t));
    fptype phi   = (thrust::get<4>(t));

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);

    // skip to our resonance function
    while(pc.funcIdx < _resonance_i)
        pc.incrementIndex();

    unsigned int pair = pc.getConstant(0);

    if(pair < 2) {
        fptype mres = pair == 0 ? m12 : m34;
        fptype d1   = pair == 0 ? m1 : m3;
        fptype d2   = pair == 0 ? m2 : m4;
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
    } else {
        fptype vecs[16];
        // TODO: What is indices[1]?
        get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
        fptype d1, d2;
        fptype mres = getmass(pair, d1, d2, vecs, m1, m2, m3, m4);
        ret         = getResonanceAmplitude(mres, d1, d2, pc);
    }

    // printf("NormLS %f, %f, %f, %f, %f \n",m12, m34, cos12, cos34, phi );
    // printf("%i, %i, %i, %i, %i \n",numLS, numSF, numAmps, offset, evtNum );
    // printf("NLS %i, %f, %f\n",_resonance_i,ret.real, ret.imag);

    // printf("m12 %f \n", m12); // %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    // printf("#Parameters %i, #LS %i, #SF %i, #AMP %i \n", indices[0], indices[3], indices[4], indices[5]);
    return ret;
}

} // namespace GooFit
