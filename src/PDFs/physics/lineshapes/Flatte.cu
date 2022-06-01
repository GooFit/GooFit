#include <goofit/PDFs/physics/lineshapes/Flatte.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include "Common.h"

namespace GooFit {

__device__ auto aSqrtTerm(const fptype &m0, const fptype &m) -> fpcomplex {
    fptype a2           = 1 - (2 * m0 / m) * (2 * m0 / m);
    fpcomplex returnVal = a2 > 0 ? fpcomplex(sqrt(a2), 0) : fpcomplex(0, sqrt(-a2));
    return returnVal;
}

__device__ auto Flatte_MINT(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    unsigned int orbital = pc.getConstant(1);
    fptype meson_radius  = pc.getConstant(2);
    fptype resmass       = pc.getParameter(0);

    fptype frFactor = 1;
    fptype rMass2   = Mpair * Mpair;

    // As far as I understand, this is only valid for the f980
    fptype gPi       = .165;
    fptype gK_by_gPi = 4.21;
    fptype gK        = gPi * gK_by_gPi;
    fptype mPi0      = .1349766;
    fptype mK0       = .497648;

    fptype mpsq = (m1 + m2) * (m1 + m2);
    fptype mmsq = (m1 - m2) * (m1 - m2);
    fptype num  = (rMass2 - mpsq) * (rMass2 - mmsq);
    // fptype num2  = (resmass*resmass - mpsq)*(resmass*resmass - mmsq);
    fptype pABSq = num / (4 * rMass2);
    // fptype prSq = fabs(num2/(4*resmass*resmass));

    fpcomplex Gpipi       = (1. / 3.) * aSqrtTerm(mPi0, Mpair) + (2. / 3.) * aSqrtTerm(mPiPlus, Mpair);
    fpcomplex GKK         = (1. / 2.) * aSqrtTerm(mK0, Mpair) + (1. / 2.) * aSqrtTerm(mKPlus, Mpair);
    fpcomplex FlatteWidth = gPi * Gpipi + gK * GKK;
    // printf("%.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Gpipi.real, Gpipi.imag, GKK.real, GKK.imag, FlatteWidth.real,
    // FlatteWidth.imag, Mpair, pABSq);

    frFactor     = BL2(pABSq * meson_radius * meson_radius, orbital);
    fpcomplex BW = sqrt(frFactor) / fpcomplex(resmass * resmass - rMass2, 0) - fpcomplex(0, 1) * resmass * FlatteWidth;

    pc.incrementIndex(1, 1, 3, 0, 1);
    return BW;
}

__device__ resonance_function_ptr ptr_to_Flatte = Flatte_MINT;

Lineshapes::Flatte::Flatte(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("Flatte", name, L, Mpair, FormFac, radius) {
    // TODO: Clean up signature
    registerParameter(mass);

    registerConstant(L);
    registerConstant(radius);

    registerFunction("ptr_to_Flatte", ptr_to_Flatte);

    initialize();
}

bool Lineshapes::Flatte::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

} // namespace GooFit
