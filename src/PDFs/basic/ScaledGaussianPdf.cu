#include "goofit/PDFs/basic/ScaledGaussianPdf.h"
#include "goofit/Variable.h"

//#include <limits>

namespace GooFit {

__device__ fptype device_ScaledGaussian(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x     = evt[0];
    fptype mean  = p[indices[1]] + p[indices[3]];
    fptype sigma = p[indices[2]] * (1 + p[indices[4]]);
    fptype ret   = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    return ret;
}

__device__ device_function_ptr ptr_to_ScaledGaussian = device_ScaledGaussian;

__host__ ScaledGaussianPdf::ScaledGaussianPdf(
    std::string n, Variable *_x, Variable *mean, Variable *sigma, Variable *delta, Variable *epsilon)
    : GooPdf(_x, n) {
    registerParameter(mean);
    registerParameter(sigma);
    registerParameter(delta);
    registerParameter(epsilon);

    std::vector<unsigned int> pindices;
    pindices.push_back(mean->getIndex());
    pindices.push_back(sigma->getIndex());
    pindices.push_back(delta->getIndex());
    pindices.push_back(epsilon->getIndex());
    GET_FUNCTION_ADDR(ptr_to_ScaledGaussian);
    initialize(pindices);
}

} // namespace GooFit
