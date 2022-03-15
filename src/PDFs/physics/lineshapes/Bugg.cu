#include <goofit/PDFs/physics/lineshapes/Bugg.h>
#include <goofit/PDFs/physics/lineshapes/Bugg3.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include "Common.h"

namespace GooFit {

__device__ auto bugg_rho2(const fptype &s, const fptype m) -> fpcomplex {
    fptype rho_squared  = 1. - 4. * m * m / s;
    fpcomplex returnVal = (rho_squared >= 0) ? fpcomplex(1, 0) : fpcomplex(0, 1);
    rho_squared         = (rho_squared >= 0) ? sqrt(rho_squared) : sqrt(-rho_squared);
    return rho_squared * returnVal;
}

__device__ auto bugg_j1(const fptype &s, const fptype m) -> fptype {
    fptype rho_pipi  = bugg_rho2(s, m).real();
    fptype returnVal = 2.;
    returnVal += (rho_pipi > 0.) ? rho_pipi * log((1. - rho_pipi) / (1. + rho_pipi)) : 0;
    return returnVal / M_PI;
}

__device__ auto bugg_Gamma_4pi(const fptype &s,
                               const fptype mpi,
                               const fptype &g_4pi,
                               const fptype &M,
                               const fptype &lambda_4pi,
                               const fptype &s0_4pi) -> fptype {
    fptype returnVal = (s < (16. * mpi * mpi)) ? 0
                                               : g_4pi * (1. / (1 + exp(lambda_4pi * (s0_4pi - s))))
                                                     / (1. / (1 + exp(lambda_4pi * (s0_4pi - M * M))));
    return returnVal;
}

// This function is an adaptation from the bugg lineshape implemented in the MINT package written by Jonas Rademacker.
// this lineshape is not tested yet!
__device__ auto bugg_MINT(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    fptype s = Mpair * Mpair;

    fptype M          = 0.953;
    fptype b1         = 1.302;
    fptype b2         = 0.340;
    fptype A          = 2.426;
    fptype g_4pi      = 0.011;
    fptype g_2K       = 0.6;
    fptype g_2eta     = 0.2;
    fptype alpha      = 1.3;
    fptype sA         = 0.41;
    fptype s0_4pi     = 7.082 / 2.845;
    fptype lambda_4pi = 2.845;

    fptype g1sq = (b1 + b2 * s) * exp(-(s - M * M) / A);
    fptype z    = bugg_j1(s, mPiPlus) - bugg_j1(M * M, mPiPlus);

    fpcomplex gamma_2pi = fpcomplex(
        g1sq * (s - sA * mPiPlus * mPiPlus) / (M * M - sA * mPiPlus * mPiPlus) * bugg_rho2(s, mPiPlus).real(), 0);
    fpcomplex gamma_2K = g_2K * g1sq * s / (M * M)
                         * exp((-1) * alpha * sqrt((s - 4. * mKPlus * mKPlus) * (s - 4. * mKPlus * mKPlus)))
                         * bugg_rho2(s, mKPlus);
    fpcomplex gamma_2eta = g_2eta * g1sq * s / (M * M)
                           * exp((-1) * alpha * sqrt((s - 4. * mEta * mEta) * (s - 4. * mEta * mEta)))
                           * bugg_rho2(s, mEta);
    fpcomplex gamma_4pi = fpcomplex(bugg_Gamma_4pi(s, mPiPlus, g_4pi, M, lambda_4pi, s0_4pi), 0);

    fpcomplex Gamma_tot = gamma_2pi + gamma_2K + gamma_2eta + gamma_4pi;

    // fpcomplex num = M * gamma_2pi; //only for elastic scattering, not production
    fpcomplex den
        = fpcomplex(M * M - s - M * g1sq * (s - sA * mPiPlus * mPiPlus) / (M * M - sA * mPiPlus * mPiPlus) * z, 0)
          - fpcomplex(0, 1) * M * Gamma_tot;
    fpcomplex returnVal = 1.0 / den;

    pc.incrementIndex(1, 0, 1, 0, 1);
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",gamma_2pi.real, gamma_2pi.imag, gamma_2K.real,
    // gamma_2K.imag, gamma_2eta.real, gamma_2eta.imag, gamma_4pi.real, gamma_4pi.imag);
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Mpair, Gamma_tot.real, Gamma_tot.imag, g1sq, z,
    // den.real, den.imag, returnVal.real, returnVal.imag);

    // the factor sqrt(1000) gives the correct result in comparison with mint2, I think its because BW/SBW
    // have a factor of sqrt(k) which these lineshapes dont have. For now this stays because it works. further
    // investigation needed.
    return returnVal * sqrt(1000.0);
}

__device__ auto bugg_MINT3(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    fptype s          = Mpair * Mpair;
    fptype M          = 0.953;
    fptype b1         = 1.302;
    fptype b2         = 0.340;
    fptype A          = 2.426;
    fptype g_4pi      = 0.011;
    fptype g_2K       = 0.6;
    fptype g_2eta     = 0.2;
    fptype alpha      = 1.3;
    fptype s0_4pi     = 7.082 / 2.845;
    fptype lambda_4pi = 2.845;
    fptype sA         = 0.41 * mPiPlus * mPiPlus;

    fptype g1sq      = M * (b1 + b2 * s) * exp(-(s - M * M) / A);
    fptype z         = bugg_j1(s, mPiPlus) - bugg_j1(M * M, mPiPlus);
    fptype adlerZero = (s - sA) / (M * M - sA);

    fptype mk4  = 4. * mKPlus * mKPlus;
    fptype me4  = 4. * mEta * mEta;
    fptype tmp1 = s > mk4 ? s - mk4 : mk4 - s;
    fptype tmp2 = s > me4 ? s - me4 : me4 - s;

    fpcomplex gamma_2pi  = fpcomplex(g1sq * adlerZero * bugg_rho2(s, mPiPlus).real(), 0);
    fpcomplex gamma_2K   = g_2K * g1sq * s / (M * M) * exp((-1) * alpha * tmp1) * bugg_rho2(s, mKPlus);
    fpcomplex gamma_2eta = g_2eta * g1sq * s / (M * M) * exp((-1) * alpha * tmp2) * bugg_rho2(s, mEta);
    fpcomplex gamma_4pi  = fpcomplex(bugg_Gamma_4pi(s, mPiPlus, g_4pi, M, lambda_4pi, s0_4pi), 0);

    fpcomplex Gamma_tot = gamma_2pi + gamma_2K + gamma_2eta + gamma_4pi;

    // fpcomplex num = M * gamma_2pi; //only for elastic scattering, not production
    fpcomplex den       = fpcomplex(M * M - s - adlerZero * g1sq * z, 0) - fpcomplex(0, 1) * Gamma_tot;
    fpcomplex returnVal = 1.0 / den;

    pc.incrementIndex(1, 0, 1, 0, 1);
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",gamma_2pi.real, gamma_2pi.imag, gamma_2K.real,
    // gamma_2K.imag, gamma_2eta.real, gamma_2eta.imag, gamma_4pi.real, gamma_4pi.imag);
    // printf("Bugg %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n",Mpair, Gamma_tot.real, Gamma_tot.imag, g1sq, z,
    // den.real, den.imag, returnVal.real, returnVal.imag);

    return returnVal;
}

__device__ resonance_function_ptr ptr_to_bugg_MINT3 = bugg_MINT3;
__device__ resonance_function_ptr ptr_to_bugg_MINT  = bugg_MINT;

Lineshapes::Bugg::Bugg(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("Bugg", name, L, Mpair, FormFac, radius) {
    // TODO: Clean up signature

    // registerConstant(radius);
    // registerConstant(L);
    // registerConstant(Mpair);
    // registerConstant(FF);

    registerFunction("ptr_to_bugg_MINT", ptr_to_bugg_MINT);

    initialize();
}

Lineshapes::Bugg3::Bugg3(
    std::string name, Variable mass, Variable width, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : Lineshape("Bugg3", name, L, Mpair, FormFac, radius) {
    // TODO: Clean up signature
    // registerConstant(radius);
    // registerConstant(L);
    // registerConstant(Mpair);
    // registerConstant(FF);

    registerFunction("ptr_to_bugg_MINT3", ptr_to_bugg_MINT3);

    initialize();
}

bool Lineshapes::Bugg::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

bool Lineshapes::Bugg3::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

} // namespace GooFit
