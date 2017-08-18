#include "goofit/PDFs/basic/BifurGaussPdf.h"

namespace GooFit {

__device__ fptype device_BifurGauss(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x          = evt[indices[2 + indices[0]]]; // why does indices recall itself?
    fptype mean       = p[indices[1]];
    fptype sigmaLeft  = p[indices[2]];
    fptype sigmaRight = p[indices[3]];

    // how to calculate the value of a bifurcated gaussian?
    fptype sigma = sigmaLeft;

    if(x > mean)
        sigma = sigmaRight;

    fptype ret = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));
    return ret;
}

__device__ device_function_ptr ptr_to_BifurGauss = device_BifurGauss;

__host__ BifurGaussPdf::BifurGaussPdf(std::string n, Variable *_x, Variable *mean, Variable *sigmaL, Variable *sigmaR)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigmaL));
    pindices.push_back(registerParameter(sigmaR));
    GET_FUNCTION_ADDR(ptr_to_BifurGauss);
    initialize(pindices);
}

// q: how shall the normalization of a bifurcated gaussian be calculated?
// a: a "sum" of two half-gaussians?
__host__ fptype BifurGaussPdf::integrate(fptype lo, fptype hi) const {
    unsigned int *indices
        = host_indices + parameters; // look at the global indexes vector starting at the parameters of this function

    fptype sL = host_params[indices[2]];
    fptype sR = host_params[indices[3]];

    fptype normL = 1. / (sqrt(2 * M_PI) * sL);
    fptype normR = 1. / (sqrt(2 * M_PI) * sR);

    return .5 * normL + .5 * normR;
}
} // namespace GooFit
