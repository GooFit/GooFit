#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/ScaledGaussianPdf.h>
#include <goofit/Variable.h>

// #include <limits>

namespace GooFit {

__device__ auto device_ScaledGaussian(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x     = RO_CACHE(evt[id]);
    fptype mean  = pc.getParameter(0) + pc.getParameter(2);
    fptype sigma = pc.getParameter(1) * (1 + pc.getParameter(3));
    fptype ret   = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    pc.incrementIndex(1, 4, 0, 1, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_ScaledGaussian = device_ScaledGaussian;

__host__ ScaledGaussianPdf::ScaledGaussianPdf(
    std::string n, Observable _x, Variable mean, Variable sigma, Variable delta, Variable epsilon)
    : GooPdf("ScaledGaussianPdf", n, _x, mean, sigma, delta, epsilon) {
    registerFunction("ptr_to_ScaledGaussian", ptr_to_ScaledGaussian);

    initialize();
}

} // namespace GooFit
