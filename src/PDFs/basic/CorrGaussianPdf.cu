#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/CorrGaussianPdf.h>

namespace GooFit {

__device__ auto device_CorrGaussian(fptype *evt, ParameterContainer &pc) -> fptype {
    int idx = pc.getObservable(0);
    int idy = pc.getObservable(1);

    fptype mean1  = pc.getParameter(0);
    fptype sigma1 = pc.getParameter(1);
    fptype mean2  = pc.getParameter(2);
    fptype sigma2 = pc.getParameter(3);
    fptype corr   = pc.getParameter(4);

    fptype x = RO_CACHE(evt[idx]);
    fptype y = RO_CACHE(evt[idy]);

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
    : GooPdf("CorrGaussianPdf", n, _x, _y, mean1, sigma1, mean2, sigma2, correlation) {
    registerFunction("ptr_to_CorrGaussian", ptr_to_CorrGaussian);

    initialize();
}

} // namespace GooFit
