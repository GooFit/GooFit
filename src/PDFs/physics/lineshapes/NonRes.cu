#include <goofit/PDFs/physics/lineshapes/NonRes.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include "Common.h"

namespace GooFit {

__device__ fpcomplex nonres_DP(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) {
    fptype meson_radius  = pc.getConstant(0);
    unsigned int orbital = pc.getConstant(1);

    fptype mumsRecoMass2 = Mpair * Mpair;

    fptype mpsq       = (m1 + m2) * (m1 + m2);
    fptype mmsq       = (m1 - m2) * (m1 - m2);
    fptype num        = (mumsRecoMass2 - mpsq) * (mumsRecoMass2 - mmsq);
    fptype pABSq      = num / (4 * mumsRecoMass2);
    fptype formfactor = sqrt(BL2(pABSq * meson_radius * meson_radius, orbital));

    pc.incrementIndex(1, 0, 2, 0, 1);
    return fpcomplex(1., 0.) * formfactor;
}

__device__ resonance_function_ptr ptr_to_NONRES_DP = nonres_DP;

Lineshapes::NonRes::NonRes(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("NonRes", name, L, Mpair, FormFac, radius) {
    // TODO: Clean up signature

    registerConstant(radius);
    registerConstant(L);

    registerFunction("ptr_to_NONRES_DP", ptr_to_NONRES_DP);

    initialize();
}

} // namespace GooFit
