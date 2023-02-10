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
    fptype s               = POW2(Mpair);
    constexpr fptype sNorm = mKPlus * mKPlus + mPiPlus * mPiPlus;
    constexpr fptype s12   = 0.23;
    constexpr fptype s32   = 0.27;

    fptype I12_adler = (s - s12) / sNorm;
    fptype I32_adler = (s - s32) / sNorm;

    fpcomplex rho1 = phsp_FOCUS(s, mKPlus, mPiPlus);
    fpcomplex rho2 = phsp_FOCUS(s, mKPlus, mEtap);
    fptype pmass   = 1.7919;
    Eigen::Array<fptype, 2, 1> coupling;
    coupling << 0.31072, -0.02323;
    fptype X = s / sNorm - 1;

    // construct KMatrix
    Eigen::Array<fptype, 2, 2> kMatrix;

    kMatrix(0, 0) = coupling(0) * coupling(0) / (pmass - s) + 0.79299 - 0.15099 * X + 0.00811 * POW2(X);
    kMatrix(1, 1) = coupling(1) * coupling(1) / (pmass - s) + 0.17054 - 0.0219 * X + 0.00085655 * POW2(X);
    kMatrix(1, 0) = coupling(1) * coupling(0) / (pmass - s) + 0.15040 - 0.038266 * X + 0.0022596 * POW2(X);
    kMatrix(0, 1) = coupling(0) * coupling(1) / (pmass - s) + 0.15040 - 0.038266 * X + 0.0022596 * POW2(X);

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

    pc.incrementIndex(1, 0, 2, 0, 1);

    if(mod == static_cast<unsigned int>(Lineshapes::FOCUS::Mod::Kpi)) {
        auto ret = (K11 - imag_i * rho2 * detK) / del;
        return ret;
    } else if(mod == static_cast<unsigned int>(Lineshapes::FOCUS::Mod::KEta)) {
        auto ret = K12 / del;
        return K12 / del;
    } else {
        return T32;
    }
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
