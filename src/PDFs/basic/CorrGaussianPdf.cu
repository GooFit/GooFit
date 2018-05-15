#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/CorrGaussianPdf.h>

namespace GooFit {

__device__ fptype device_CorrGaussian(fptype *evt, ParameterContainer &pc) {
    int idx = pc.getObservable(0);
    int idy = pc.getObservable(1);

    fptype mean1  = pc.getParameter(0);
    fptype sigma1 = pc.getParameter(1);
    fptype mean2  = pc.getParameter(2);
    fptype sigma2 = pc.getParameter(3);
    fptype corr   = pc.getParameter(4);

    fptype x = evt[idx];
    fptype y = evt[idy];

    fptype x_dist = (x - mean1) / sigma1;
    sigma2 *= (1.0 + corr * x_dist * x_dist);
    fptype y_dist = (y - mean2) * (sigma2 == 0.0 ? 0.0 : (1.0 / sigma2));
    fptype ret    = exp(-0.5 * (x_dist * x_dist + y_dist * y_dist));

    // if ((gpuDebug & 1) && (THREADIDX == 60)) printf("CorrGauss: %f %f %f %f %f %f %f %f %f %.12f\n", x, y, mean1,
    // sigma1, mean2, sigma2, corr, x_dist, y_dist, ret);
    // if ((gpuDebug & 1) && (THREADIDX == 60)) printf("[%i, %i] [%i, %i] CorrGauss: %f %f %f\n", BLOCKIDX, THREADIDX,
    // gridDim.x, BLOCKDIM, x, y, ret);
    // printf("CorrGauss: %i %i %i %f %f %f %f\n", indices[2 + indices[0]], indices[3 + indices[0]], indices[0], x, y,
    // mean1, mean2);

    pc.incrementIndex(1, 5, 0, 2, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_CorrGaussian = device_CorrGaussian;

__host__ CorrGaussianPdf::CorrGaussianPdf(std::string n,
                                          Observable _x,
                                          Observable _y,
                                          Variable mean1,
                                          Variable sigma1,
                                          Variable mean2,
                                          Variable sigma2,
                                          Variable correlation)
    : GooPdf(n, _x, _y) {
    // registerObservable(_x);
    // registerObservable(_y);

    registerParameter(mean1);
    registerParameter(sigma1);
    registerParameter(mean2);
    registerParameter(sigma2);
    registerParameter(correlation);

    initialize();
}

__host__ void CorrGaussianPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_CorrGaussian);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_CorrGaussian");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

} // namespace GooFit
