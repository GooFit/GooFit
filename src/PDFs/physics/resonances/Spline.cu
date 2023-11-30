#include <goofit/PDFs/physics/resonances/Spline.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/detail/ComplexUtils.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ fptype cDeriatives[2 * 100];

__device__ auto cubicspline(fptype m13, fptype m23, fptype m12, ParameterContainer &pc) -> fpcomplex {
    fpcomplex ret(0, 0);
    unsigned int cyclic_index        = pc.getConstant(0);
    unsigned int doSwap              = pc.getConstant(1);
    const unsigned int nKnobs        = pc.getConstant(2);
    const unsigned int linear        = pc.getConstant(3);
    unsigned int swave_const_idx     = 4; // num consts before swave constants
    unsigned int idx                 = 4; // Next index
    unsigned int i                   = 0;
    idx += 2 * nKnobs;
    fptype mAC = m13, mBC = m23;
    switch(cyclic_index) {
    case PAIR_13:
        mAC = m13;
        mBC = m23;
        break;
    case PAIR_23:
        mAC = m23;
        mBC = m12;
        break;
    }

    int khiAC = 0, khiBC = 0;
    fptype dmKK, aa, bb, aa3, bb3;
    unsigned int timestorun = 1 + doSwap;

    // Run 0 and/or 1
    for(i = 0; i < timestorun; i++) {
        // Find the knots we are between
        while(khiAC < nKnobs) {
            if(mAC < pc.getConstant(swave_const_idx + khiAC))
                break;
            khiAC++;
        }

        // Quit this iteration if outside
        if(khiAC > 0 && khiAC < nKnobs) {
            unsigned int kloAC          = khiAC - 1; //, kloAC = khiAC -1;
            unsigned int twokloAC       = kloAC + kloAC;
            unsigned int twokhiAC       = khiAC + khiAC;
            fptype pwa_coefs_real_kloAC = pc.getParameter(twokloAC);
            fptype pwa_coefs_real_khiAC = pc.getParameter(twokhiAC);
            fptype pwa_coefs_imag_kloAC = pc.getParameter(twokloAC + 1);
            fptype pwa_coefs_imag_khiAC = pc.getParameter(twokhiAC + 1);

            fptype pwa_coefs_prime_real_kloAC = cDeriatives[twokloAC];
            fptype pwa_coefs_prime_real_khiAC = cDeriatives[twokhiAC];
            fptype pwa_coefs_prime_imag_kloAC = cDeriatives[twokloAC + 1];
            fptype pwa_coefs_prime_imag_khiAC = cDeriatives[twokhiAC + 1];

            dmKK = pc.getConstant(swave_const_idx + khiAC) - pc.getConstant(swave_const_idx + kloAC);
            aa   = (pc.getConstant(swave_const_idx+khiAC) - mAC) / dmKK;
            bb   = 1 - aa;
            aa3  = aa * aa * aa;
            bb3  = bb * bb * bb;
            
            if(linear){
                ret.real(ret.real() + aa * pwa_coefs_real_kloAC + bb * pwa_coefs_real_khiAC);
                ret.imag(ret.imag() + aa * pwa_coefs_imag_kloAC + bb * pwa_coefs_imag_khiAC);
            }else{
                 ret.real(ret.real() + aa * pwa_coefs_real_kloAC + bb * pwa_coefs_real_khiAC
                     + ((aa3 - aa) * pwa_coefs_prime_real_kloAC + (bb3 - bb) * pwa_coefs_prime_real_khiAC)
                           * (dmKK * dmKK) / 6.0);
                 ret.imag(ret.imag() + aa * pwa_coefs_imag_kloAC + bb * pwa_coefs_imag_khiAC
                     + ((aa3 - aa) * pwa_coefs_prime_imag_kloAC + (bb3 - bb) * pwa_coefs_prime_imag_khiAC)
                           * (dmKK * dmKK) / 6.0);
            }
        }

        khiAC = khiBC;
        mAC   = mBC;
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
               bool symmDP,
                bool linear)
    : ResonancePdf("Spline", name, ar, ai) {
    const unsigned int nKnobs = HH_bin_limits.size();

    registerConstant(cyc);
    registerConstant(symmDP);
    registerConstant(nKnobs);
    registerConstant(linear);

    for(int i = 0; i < pwa_coefs_reals.size(); i++) {
        registerConstant(HH_bin_limits[i]);
        registerParameter(pwa_coefs_reals[i]);
        registerParameter(pwa_coefs_imags[i]);
    }

    registerFunction("ptr_to_SPLINE", ptr_to_SPLINE);

    recalculateCache();
}

__host__ void Spline::recalculateCache() const {
    unsigned int swave_const_idx     = 4; // num consts before swave constants
    auto params           = getParameters();
    const unsigned nKnobs = constantsList[2];
    std::vector<fptype> x(nKnobs);
    std::vector<fpcomplex> y(nKnobs);
    unsigned int i = 0;
    for(auto v = params.begin(); v != params.end(); ++v, ++i) {
        unsigned int idx = i / 2;
        fptype value     = parametersList[i];
        if(i % 2 == 0) {
            x[idx] = constantsList[swave_const_idx + idx];
            y[idx].real(value);
        } else
            y[idx].imag(value);
    }

    std::vector<fptype> y2_flat = flatten(complex_derivative(x, y));

    MEMCPY_TO_SYMBOL(cDeriatives, y2_flat.data(), 2 * nKnobs * sizeof(fptype), 0, cudaMemcpyHostToDevice);
}

} // namespace Resonances
} // namespace GooFit
