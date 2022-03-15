#include <goofit/PDFs/physics/lineshapes/LASS.h>
#include <goofit/PDFs/physics/lineshapes/RBW.h>
#include <goofit/PDFs/physics/lineshapes/SBW.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include "Common.h"

namespace GooFit {

// This function is modeled after BW_BW::getVal() in BW_BW.cpp from the MINT package written by Jonas Rademacker.
__device__ auto BW(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    fptype resmass       = pc.getParameter(0);
    fptype reswidth      = pc.getParameter(1);
    unsigned int orbital = pc.getConstant(1);
    unsigned int FF      = pc.getConstant(2);
    fptype meson_radius  = pc.getConstant(3);

    const unsigned int to2Lplus1 = 2 * orbital + 1;

    fptype mass          = resmass;
    fptype width         = reswidth;
    fptype mumsRecoMass2 = Mpair * Mpair;

    fptype mpsq        = (m1 + m2) * (m1 + m2);
    fptype mmsq        = (m1 - m2) * (m1 - m2);
    fptype num         = (mumsRecoMass2 - mpsq) * (mumsRecoMass2 - mmsq);
    fptype num2        = (mass * mass - mpsq) * (mass * mass - mmsq);
    fptype pABSq       = num / (4 * mumsRecoMass2);
    fptype prSqForGofM = num2 / (4 * mass * mass);
    fptype prSq2       = prSqForGofM < 0 ? 0 : prSqForGofM;
    prSqForGofM        = fabs(prSqForGofM);

    fptype pratio = sqrt(pABSq / prSqForGofM);

    fptype pratio_to_2Jplus1 = 1;

    for(int i = 0; i < to2Lplus1; i++) {
        pratio_to_2Jplus1 *= pratio;
    }

    fptype mratio   = mass / Mpair;
    fptype r2       = meson_radius * meson_radius;
    fptype thisFR   = BL_PRIME(pABSq * r2, prSqForGofM * r2, orbital);
    fptype frFactor = 1;

    if(0 != orbital && 0 != FF) {
        frFactor = (FF == 1 ? BL(pABSq * r2, orbital) : BL_PRIME(pABSq * r2, prSq2 * r2, orbital));
        frFactor = (FF == 3 ? BL2(pABSq * r2, orbital) : frFactor);
    }

    fptype GofM = width * pratio_to_2Jplus1 * mratio * thisFR;

    fptype gamma = mass * sqrt((mass * mass + width * width));
    fptype k     = (2.0 * sqrt(2.0) / M_PI) * mass * width * gamma
               / sqrt(mass * mass + gamma); // Note added additional factor of 2*sqrt(2)/PI here so results are
                                            // comparable to MINT3. MINT2 doesn't have include this.

    fpcomplex BW(mass * mass - mumsRecoMass2, mass * GofM);
    fptype den = (mass * mass - mumsRecoMass2) * (mass * mass - mumsRecoMass2) + mass * GofM * mass * GofM;

    fpcomplex ret = (sqrt(k * frFactor)) / den * BW;

    pc.incrementIndex(1, 2, 4, 0, 1);

    return ret;
}

// This function is modeled after SBW from the MINT package written by Jonas Rademacker.
__device__ auto SBW(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    fptype resmass       = pc.getParameter(0);
    fptype reswidth      = pc.getParameter(1);
    unsigned int orbital = pc.getConstant(1);
    unsigned int FF      = pc.getConstant(2);
    fptype meson_radius  = pc.getConstant(3);

    fptype mass          = resmass;
    fptype width         = reswidth;
    fptype mumsRecoMass2 = Mpair * Mpair;

    fptype mpsq  = (m1 + m2) * (m1 + m2);
    fptype mmsq  = (m1 - m2) * (m1 - m2);
    fptype num   = (mumsRecoMass2 - mpsq) * (mumsRecoMass2 - mmsq);
    fptype num2  = (mass * mass - mpsq) * (mass * mass - mmsq);
    fptype pABSq = num / (4 * mumsRecoMass2);
    fptype prSq  = num2 / (4 * mass * mass);
    fptype prSq2 = prSq < 0 ? 0 : prSq;
    prSq         = fabs(prSq);

    fptype r2       = meson_radius * meson_radius;
    fptype frFactor = 1;

    if(0 != orbital && 0 != FF) {
        frFactor = (FF == 1 ? BL(pABSq * r2, orbital) : BL_PRIME(pABSq * r2, prSq2 * r2, orbital));
        frFactor = (FF == 3 ? BL2(pABSq * r2, orbital) : frFactor);
    }

    fptype GofM = width;

    fptype gamma = sqrt(mass * mass * (mass * mass + width * width));
    fptype k     = mass * width * gamma / sqrt(mass * mass + gamma);

    fpcomplex BW(mass * mass - mumsRecoMass2, mass * GofM);
    fptype den = (mass * mass - mumsRecoMass2) * (mass * mass - mumsRecoMass2) + mass * GofM * mass * GofM;

    fpcomplex ret = (sqrt(k * frFactor)) / den * BW;

    pc.incrementIndex(1, 2, 4, 0, 1);

    // printf("m1, m2, Mpair, GofM, pABSq , prSq, FF, ret.real, ret.imag\n");
    // printf("SBW %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g\n", m1, m2, Mpair, GofM, pABSq, prSq, frFactor,
    // ret.real, ret.imag );
    return ret;
}

__device__ auto lass_MINT(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);
    fptype rMass2   = Mpair * Mpair;

    fptype a = 2.07;
    fptype r = 3.32;

    fptype mpsq  = (m1 + m2) * (m1 + m2);
    fptype mmsq  = (m1 - m2) * (m1 - m2);
    fptype num   = (rMass2 - mpsq) * (rMass2 - mmsq);
    fptype num2  = (resmass * resmass - mpsq) * (resmass * resmass - mmsq);
    fptype pABSq = num / (4 * rMass2);
    fptype prSq  = fabs(num2 / (4 * resmass * resmass));

    fptype y          = 2.0 * a * sqrt(pABSq);
    fptype x          = 2.0 + a * r * pABSq;
    fptype cotDeltaBg = x / y;
    fpcomplex phaseshift((cotDeltaBg * cotDeltaBg - 1) / (1 + cotDeltaBg * cotDeltaBg),
                         2 * cotDeltaBg / (1 + cotDeltaBg * cotDeltaBg));
    // (cotDeltaBg*cotDeltaBg-1)/(1+cotDeltaBg*cotDeltaBg) = cos(2*delta)     2*cotDeltaBg / ( 1 +
    // cotDeltaBg*cotDeltaBg) = sin(2*delta)
    fpcomplex den(sqrt(pABSq) * cotDeltaBg, (-1.) * sqrt(pABSq));
    fptype SF           = Mpair * sqrt(prSq) / (resmass * resmass * reswidth);
    fpcomplex BG        = SF / den;
    fpcomplex returnVal = BG + phaseshift * BW(Mpair, m1, m2, pc);

    pc.incrementIndex(1, 2, 4, 0, 1);
    // The call to BW does the increment for us

    return returnVal;
}

__device__ resonance_function_ptr ptr_to_BW_DP4 = BW;
__device__ resonance_function_ptr ptr_to_SBW    = SBW;
__device__ resonance_function_ptr ptr_to_lass   = lass_MINT;

Lineshapes::RBW::RBW(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("RBW", name, L, Mpair, FormFac, radius) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(L);
    registerConstant(enum_to_underlying(FormFac));
    registerConstant(radius);

    registerFunction("ptr_to_BW_DP4", ptr_to_BW_DP4);

    initialize();
}

bool Lineshapes::RBW::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

Lineshapes::SBW::SBW(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("SBW", name, L, Mpair, FormFac, radius) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(L);
    registerConstant(enum_to_underlying(FormFac));
    registerConstant(radius);

    registerFunction("ptr_to_SBW", ptr_to_SBW);

    initialize();
}

bool Lineshapes::SBW::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

Lineshapes::LASS::LASS(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("LASS", name, L, Mpair, FormFac, radius) {
    // This one must match RBW exactly, because it calls the same function at one point

    registerParameter(mass);
    registerParameter(width);

    registerConstant(L);
    registerConstant(enum_to_underlying(FormFac));
    registerConstant(radius);

    registerFunction("ptr_to_lass", ptr_to_lass);

    initialize();
}

bool Lineshapes::LASS::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

} // namespace GooFit
