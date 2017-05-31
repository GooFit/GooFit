#include "goofit/PDFs/basic/ExpGausPdf.h"

namespace GooFit {

__device__ fptype device_ExpGaus(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x     = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])];
    fptype mean  = RO_CACHE(p[RO_CACHE(indices[1])]);
    fptype sigma = RO_CACHE(p[RO_CACHE(indices[2])]);
    fptype alpha = RO_CACHE(p[RO_CACHE(indices[3])]);

    fptype ret    = 0.5 * alpha;
    fptype exparg = ret * (2 * mean + alpha * sigma * sigma - 2 * x);
    fptype erfarg = (mean + alpha * sigma * sigma - x) / (sigma * 1.4142135623);

    ret *= exp(exparg);
    ret *= erfc(erfarg);

    return ret;
}

__device__ device_function_ptr ptr_to_ExpGaus = device_ExpGaus;

ExpGausPdf::ExpGausPdf(std::string n, Variable *_x, Variable *mean, Variable *sigma, Variable *tau)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigma));
    pindices.push_back(registerParameter(tau));
    GET_FUNCTION_ADDR(ptr_to_ExpGaus);
    initialize(pindices);
}

} // namespace GooFit
