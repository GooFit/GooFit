#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/BifurGaussPdf.h>

namespace GooFit {

__device__ auto device_BifurGauss(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x          = RO_CACHE(evt[id]);
    fptype mean       = pc.getParameter(0);
    fptype sigmaLeft  = pc.getParameter(1);
    fptype sigmaRight = pc.getParameter(2);

    // how to calculate the value of a bifurcated gaussian?
    fptype sigma = sigmaLeft;

    if(x > mean)
        sigma = sigmaRight;

    pc.incrementIndex(1, 3, 0, 1, 1);

    fptype ret = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));
    return ret;
}

__device__ device_function_ptr ptr_to_BifurGauss = device_BifurGauss;

__host__ BifurGaussPdf::BifurGaussPdf(std::string n, Observable _x, Variable mean, Variable sigmaL, Variable sigmaR)
    : GooPdf("BifurGaussPdf", n, _x, mean, sigmaL, sigmaR) {
    registerFunction("ptr_to_BifurGauss", ptr_to_BifurGauss);

    initialize();
}

// q: how shall the normalization of a bifurcated gaussian be calculated?
// a: a "sum" of two half-gaussians?
__host__ auto BifurGaussPdf::integrate(fptype lo, fptype hi) const -> fptype {
    fptype sL = host_parameters[parametersIdx + 2];
    fptype sR = host_parameters[parametersIdx + 3];

    fptype normL = 1. / (sqrt(2 * M_PI) * sL);
    fptype normR = 1. / (sqrt(2 * M_PI) * sR);

    return .5 * normL + .5 * normR;
}
} // namespace GooFit
