#include "goofit/PDFs/basic/CorrGaussianPdf.h"

namespace GooFit {

__device__ fptype device_CorrGaussian(fptype *evt, ParameterContainer &pc) {
    int idx = pc.observables[pc.observableIdx + 1];
    int idy = pc.observables[pc.observableIdx + 2];

    fptype mean1  = pc.parameters[pc.parameterIdx + 1];
    fptype sigma1 = pc.parameters[pc.parameterIdx + 2];
    fptype mean2  = pc.parameters[pc.parameterIdx + 3];
    fptype sigma2 = pc.parameters[pc.parameterIdx + 4];
    fptype corr   = pc.parameters[pc.parameterIdx + 5];

    fptype x      = evt[idx];
    fptype y      = evt[idy];

    fptype x_dist = (x - mean1) / sigma1;
    sigma2 *= (1 + corr * x_dist * x_dist);
    fptype y_dist = (y - mean2) * (sigma2 == 0 ? 0 : (1.0 / sigma2));
    fptype ret    = exp(-0.5 * (x_dist * x_dist + y_dist * y_dist));

    // if ((gpuDebug & 1) && (THREADIDX == 60)) printf("CorrGauss: %f %f %f %f %f %f %f %f %f %.12f\n", x, y, mean1,
    // sigma1, mean2, sigma2, corr, x_dist, y_dist, ret);
    // if ((gpuDebug & 1) && (THREADIDX == 60)) printf("[%i, %i] [%i, %i] CorrGauss: %f %f %f\n", BLOCKIDX, THREADIDX,
    // gridDim.x, BLOCKDIM, x, y, ret);
    // printf("CorrGauss: %i %i %i %f %f %f %f\n", indices[2 + indices[0]], indices[3 + indices[0]], indices[0], x, y,
    // mean1, mean2);

    pc.incrementIndex (1, 5, 0, 2, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_CorrGaussian = device_CorrGaussian;

__host__ CorrGaussianPdf::CorrGaussianPdf(std::string n,
                                          Variable *_x,
                                          Variable *_y,
                                          Variable *mean1,
                                          Variable *sigma1,
                                          Variable *mean2,
                                          Variable *sigma2,
                                          Variable *correlation)
    : GooPdf(_x, n) {
    // registerObservable(_x);
    registerObservable(_y);

    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(mean1));
    pindices.push_back(registerParameter(sigma1));
    pindices.push_back(registerParameter(mean2));
    pindices.push_back(registerParameter(sigma2));
    pindices.push_back(registerParameter(correlation));

    GET_FUNCTION_ADDR(ptr_to_CorrGaussian);
    initialize(pindices);
}

__host__ void CorrGaussianPdf::recursiveSetIndices () {
    GET_FUNCTION_ADDR(ptr_to_CorrGaussian);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_CorrGaussian");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;

    populateArrays ();
}

} // namespace GooFit
