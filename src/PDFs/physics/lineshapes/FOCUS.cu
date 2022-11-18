#include <goofit/PDFs/physics/lineshapes/FOCUS.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "Common.h"

namespace GooFit {

__device__ auto phsp_FOCUS(fptype s, fptype m0, fptype m1) -> fpcomplex {
    fptype mp    = (m0 + m1);
    fptype mm    = (m0 - m1);
    fpcomplex a2 = (1.0 - mp * mp / s) * (1.0 - mm * mm / s);
    return sqrt(a2);
}

__device__ auto FOCUSFunction(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    unsigned int mod = pc.getConstant(1);

    fpcomplex imag_i(0, 1); // imaginary i
    fptype s = POW2(Mpair);
    // printf("FOCUS func args, Mpair:%.7g, m1:%.7g, m2:%.7g, s:%.7g\n", Mpair, m1, m2, s);
    //  mKPlus, mPiPlus, mEtap
    constexpr fptype sNorm = mKPlus * mKPlus + mPiPlus * mPiPlus;
    constexpr fptype s12   = 0.23;
    constexpr fptype s32   = .27;

    fptype I12_adler = (s - s12) / sNorm;
    fptype I32_adler = (s - s32) / sNorm;

    fpcomplex rho1 = phsp_FOCUS(s, mKPlus, mPiPlus);
    fpcomplex rho2 = phsp_FOCUS(s, mKPlus, mEtap);
    // printf("rho1 re: %.7g,rho1 im: %.7g , rho2 re:%.7g, rho2
    // im:%.7g\n",rho1.real(),rho1.imag(),rho2.real(),rho2.imag());
    fptype pmass = 1.7919;
    Eigen::Array<fptype, 2, 1> coupling;
    coupling << 0.31072, -0.02323;
    fptype X = s / sNorm - 1;

    // constructKMatrix
    Eigen::Array<fptype, 2, 2> kMatrix;

    kMatrix(0, 0) = coupling(0) * coupling(0) / (pmass - s) + 0.79299 - 0.15099 * X + 0.00811 * POW2(X);
    kMatrix(1, 1) = coupling(1) * coupling(1) / (pmass - s) + 0.17054 - 0.0219 * X + 0.00085655 * POW2(X);
    kMatrix(1, 0) = coupling(1) * coupling(0) / (pmass - s) + 0.15040 - 0.038266 * X + 0.0022596 * POW2(X);
    kMatrix(0, 1) = coupling(0) * coupling(1) / (pmass - s) + 0.15040 - 0.038266 * X + 0.0022596 * POW2(X);

    // printf("FOCUS kMatrix elements, (0,0):%.7g, (1,1):%.7g, (1,0):%.7g,
    // (0,1):%.7g\n",kMatrix(0,0),kMatrix(1,1),kMatrix(1,0),kMatrix(0,1));

    fptype K11 = I12_adler * kMatrix(0, 0);
    fptype K12 = I12_adler * kMatrix(0, 1);
    fptype K22 = I12_adler * kMatrix(1, 1);

    fptype K32 = I32_adler * (-0.22147 + 0.026637 * X - 0.00092057 * POW2(X));

    fptype detK   = K11 * K22 - K12 * K12;
    fpcomplex del = 1. - rho1 * rho2 * detK - imag_i * (rho1 * K11 + rho2 * K22);

    fpcomplex T11 = 1. - imag_i * rho2 * K22;
    fpcomplex T22 = 1. - imag_i * rho1 * K11;
    fpcomplex T12 = imag_i * rho2 * K12;

    fpcomplex T32 = 1. / (1. - imag_i * rho2 * K32);

    // printf("K11:%.7g, K12:%.7g, K22:%.7g, K32:%.7g, detK:%.7g, T11 re:%.7g, T11 im:%.7g, T22 re:%.7g, T22 im:%.7g,
    // T12 re:%.7g, T12 im:%.7g, T32 re:%.7g, T32 im:%.7g, del re:%.7g, del im:%.7g\n",
    // K11,K12,K22,K32,detK,T11.real(),T11.imag(),T22.real(),
    // T22.imag(),T12.real(),T12.imag(),T32.real(),T32.imag(),del.real(), del.imag());

    pc.incrementIndex(1, 0, 2, 0, 1);

    if(mod == static_cast<unsigned int>(Lineshapes::FOCUS::Mod::Kpi)) {
        auto ret = (K11 - imag_i * rho2 * detK) / del;
        // printf("FOCUS Kpi modification return value real:%.7g, imag:%.7g\n",ret.real(),ret.imag());
        return ret;
    } else if(mod == static_cast<unsigned int>(Lineshapes::FOCUS::Mod::KEta)) {
        auto ret = K12 / del;
        // printf("FOCUS KEta modification return value real:%.7g, imag:%.7g\n",ret.real(),ret.imag());
        return K12 / del;
    } else { /*if(mod==Lineshapes::FOCUS::Mod::I32)*/

        // printf("FOCUS I32 return value real:%.7g, imag:%.7g\n",T32.real(),T32.imag());
        return T32;
    }
    // printf("Didn't modify anything. Returning zero.");
    return {0., 0.};
}

__device__ resonance_function_ptr ptr_to_FOCUS = FOCUSFunction;

Lineshapes::FOCUS::FOCUS(std::string name,
                         Mod mod,
                         Variable mass,
                         Variable width,
                         unsigned int L,
                         unsigned int Mpair,
                         FF FormFac,
                         fptype radius)
    : Lineshape("FOCUS", name, L, Mpair, FormFac, radius) {
    // TODO: Clean up signature
    registerConstant(static_cast<unsigned int>(mod));

    registerFunction("ptr_to_FOCUS", ptr_to_FOCUS);

    initialize();
}

bool Lineshapes::FOCUS::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

} // namespace GooFit
