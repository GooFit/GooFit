#include <goofit/PDFs/physics/resonances/Spline.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/detail/ComplexUtils.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ fptype cDeriatives[2 * 100];

__device__ auto cubicspline(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    fpcomplex ret(0, 0);
    unsigned int cyclic_index        = pc.getConstant(0);
    unsigned int doSwap              = pc.getConstant(1);
    const unsigned int nKnobs        = pc.getConstant(2);
    unsigned int idx                 = 3; // Next index
    unsigned int i                   = 0;
    const unsigned int pwa_coefs_idx = idx;
    idx += 2 * nKnobs;
    fptype mAB = m12, mAC = m13;
    switch(cyclic_index) {
    case PAIR_13:
        mAB = m13;
        mAC = m12;
        break;
    case PAIR_23:
        mAB = m23;
        mAC = m12;
        break;
    }

    int khiAB = 0, khiAC = 0;
    fptype dmKK, aa, bb, aa3, bb3;
    unsigned int timestorun = 1 + doSwap;

    // Run 0 and/or 1
    for(i = 0; i < timestorun; i++) {
        // Find the knots we are between
        while(khiAB < nKnobs) {
            if(mAB < pc.getConstant(3 + khiAB))
                break;
            khiAB++;
        }

        // Quit this iteration if outside
        if(khiAB > 0 && khiAB < nKnobs) {
            unsigned int kloAB          = khiAB - 1; //, kloAC = khiAC -1;
            unsigned int twokloAB       = kloAB + kloAB;
            unsigned int twokhiAB       = khiAB + khiAB;
            fptype pwa_coefs_real_kloAB = pc.getParameter(twokloAB);
            fptype pwa_coefs_real_khiAB = pc.getParameter(twokhiAB);
            fptype pwa_coefs_imag_kloAB = pc.getParameter(twokloAB + 1);
            fptype pwa_coefs_imag_khiAB = pc.getParameter(twokhiAB + 1);

            fptype pwa_coefs_prime_real_kloAB = cDeriatives[twokloAB];
            fptype pwa_coefs_prime_real_khiAB = cDeriatives[twokhiAB];
            fptype pwa_coefs_prime_imag_kloAB = cDeriatives[twokloAB + 1];
            fptype pwa_coefs_prime_imag_khiAB = cDeriatives[twokhiAB + 1];

            dmKK = pc.getConstant(3 + khiAB) - pc.getConstant(3 + kloAB);
            aa   = (pc.getConstant(3 + khiAB) - mAB) / dmKK;
            bb   = 1 - aa;
            aa3  = aa * aa * aa;
            bb3  = bb * bb * bb;
            ret.real(ret.real() + aa * pwa_coefs_real_kloAB + bb * pwa_coefs_real_khiAB
                     + ((aa3 - aa) * pwa_coefs_prime_real_kloAB + (bb3 - bb) * pwa_coefs_prime_real_khiAB)
                           * (dmKK * dmKK) / 6.0);
            ret.imag(ret.imag() + aa * pwa_coefs_imag_kloAB + bb * pwa_coefs_imag_khiAB
                     + ((aa3 - aa) * pwa_coefs_prime_imag_kloAB + (bb3 - bb) * pwa_coefs_prime_imag_khiAB)
                           * (dmKK * dmKK) / 6.0);
        }

        khiAB = khiAC;
        mAB   = mAC;
    }
    return ret;
}

__device__ resonance_function_ptr ptr_to_SPLINE = cubicspline;

namespace Resonances {

Spline::Spline(std::string name,
               Variable ar,
               Variable ai,
               std::vector<fptype> &HH_bin_limits,
               std::vector<Variable> &pwa_coefs_reals,
               std::vector<Variable> &pwa_coefs_imags,
               unsigned int cyc,
               bool symmDP)
    : ResonancePdf("Spline", name, ar, ai) {
    const unsigned int nKnobs = HH_bin_limits.size();

    registerConstant(cyc);
    registerConstant(symmDP);
    registerConstant(nKnobs);

    for(int i = 0; i < pwa_coefs_reals.size(); i++) {
        registerConstant(HH_bin_limits[i]);
        registerParameter(pwa_coefs_reals[i]);
        registerParameter(pwa_coefs_imags[i]);
    }

    registerFunction("ptr_to_SPLINE", ptr_to_SPLINE);

    recalculateCache();
}

__host__ void Spline::recalculateCache() const {
    auto params           = getParameters();
    const unsigned nKnobs = constantsList[2];
    std::vector<fptype> x(nKnobs);
    std::vector<fpcomplex> y(nKnobs);
    unsigned int i = 0;
    for(auto v = params.begin(); v != params.end(); ++v, ++i) {
        unsigned int idx = i / 2;
        fptype value     = parametersList[i];
        if(i % 2 == 0) {
            x[idx] = constantsList[3 + idx];
            y[idx].real(value);
        } else
            y[idx].imag(value);
    }

    std::vector<fptype> y2_flat = flatten(complex_derivative(x, y));

    MEMCPY_TO_SYMBOL(cDeriatives, y2_flat.data(), 2 * nKnobs * sizeof(fptype), 0, cudaMemcpyHostToDevice);
}

} // namespace Resonances
} // namespace GooFit
