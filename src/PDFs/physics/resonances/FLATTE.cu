#include <goofit/PDFs/physics/resonances/FLATTE.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ void
getAmplitudeCoefficients(fpcomplex a1, fpcomplex a2, fptype &a1sq, fptype &a2sq, fptype &a1a2real, fptype &a1a2imag) {
    // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
    a1sq = thrust::norm(a1);
    a2sq = thrust::norm(a2);
    a1 *= conj(a2);
    a1a2real = a1.real();
    a1a2imag = a1.imag();
}

__device__ auto flatte(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    // indices[1] is unused constant index, for consistency with other function types.
    fptype resmass            = pc.getParameter(0);
    fptype g1                 = pc.getParameter(1);
    fptype g2                 = pc.getParameter(2) * g1;
    unsigned int cyclic_index = pc.getConstant(0);
    unsigned int doSwap       = pc.getConstant(1);

    fptype pipmass = 0.13957018;
    fptype pi0mass = 0.1349766;
    fptype kpmass  = 0.493677;
    fptype k0mass  = 0.497614;

    fptype twopimasssq  = 4 * pipmass * pipmass;
    fptype twopi0masssq = 4 * pi0mass * pi0mass;
    fptype twokmasssq   = 4 * kpmass * kpmass;
    fptype twok0masssq  = 4 * k0mass * k0mass;

    fpcomplex ret(0., 0.);
    for(int i = 0; i < 1 + doSwap; i++) {
        fptype rhopipi_real = 0, rhopipi_imag = 0;
        fptype rhokk_real = 0, rhokk_imag = 0;

        fptype s = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

        if(s >= twopimasssq)
            rhopipi_real += (2. / 3) * sqrt(1 - twopimasssq / s); // Above pi+pi- threshold
        else
            rhopipi_imag += (2. / 3) * sqrt(-1 + twopimasssq / s);
        if(s >= twopi0masssq)
            rhopipi_real += (1. / 3) * sqrt(1 - twopi0masssq / s); // Above pi0pi0 threshold
        else
            rhopipi_imag += (1. / 3) * sqrt(-1 + twopi0masssq / s);
        if(s >= twokmasssq)
            rhokk_real += 0.5 * sqrt(1 - twokmasssq / s); // Above K+K- threshold
        else
            rhokk_imag += 0.5 * sqrt(-1 + twokmasssq / s);
        if(s >= twok0masssq)
            rhokk_real += 0.5 * sqrt(1 - twok0masssq / s); // Above K0K0 threshold
        else
            rhokk_imag += 0.5 * sqrt(-1 + twok0masssq / s);
        fptype A = (resmass * resmass - s) + resmass * (rhopipi_imag * g1 + rhokk_imag * g2);
        fptype B = resmass * (rhopipi_real * g1 + rhokk_real * g2);
        fptype C = 1.0 / (A * A + B * B);
        fpcomplex retur(A * C, B * C);
        ret += retur;
        if(doSwap) {
            fptype swpmass = m12;
            m12            = m13;
            m13            = swpmass;
        }
    }

    pc.incrementIndex(1, 3, 2, 0, 1);
    return ret;
}

__device__ resonance_function_ptr ptr_to_FLATTE = flatte;

namespace Resonances {

FLATTE::FLATTE(std::string name,
               Variable ar,
               Variable ai,
               Variable mean,
               Variable g1,
               Variable rg2og1,
               unsigned int cyc,
               bool symmDP)
    : ResonancePdf("FLATTE", name, ar, ai) {
    registerParameter(mean);
    registerParameter(g1);
    registerParameter(rg2og1);

    registerConstant(cyc);
    registerConstant(symmDP);

    registerFunction("ptr_to_FLATTE", ptr_to_FLATTE);
}

} // namespace Resonances
} // namespace GooFit
