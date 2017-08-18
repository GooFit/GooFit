#include "goofit/PDFs/basic/GaussianPdf.h"

namespace GooFit {

__device__ fptype device_Gaussian(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x     = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])];
    fptype mean  = RO_CACHE(p[RO_CACHE(indices[1])]);
    fptype sigma = RO_CACHE(p[RO_CACHE(indices[2])]);

    fptype ret = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    return ret;
}

__device__ device_function_ptr ptr_to_Gaussian = device_Gaussian;

__host__ GaussianPdf::GaussianPdf(std::string n, Variable *_x, Variable *mean, Variable *sigma)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigma));
    GET_FUNCTION_ADDR(ptr_to_Gaussian);
    initialize(pindices);
}

__host__ fptype GaussianPdf::integrate(fptype lo, fptype hi) const {
    static const fptype rootPi = sqrt(atan2(0.0, -1.0));

    unsigned int *indices = host_indices + parameters;

    // Integral over all R.
    fptype sigma = host_params[indices[2]];
    sigma *= root2 * rootPi;
    return sigma;
}

} // namespace GooFit
