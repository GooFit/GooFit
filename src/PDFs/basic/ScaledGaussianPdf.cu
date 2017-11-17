#include <goofit/PDFs/basic/ScaledGaussianPdf.h>
#include <goofit/Variable.h>

//#include <limits>

namespace GooFit {

__device__ fptype device_ScaledGaussian(fptype *evt, ParameterContainer &pc) {
    int id = pc.getObservable(0);

    fptype x     = evt[id];
    fptype mean  = pc.getParameter(0) + pc.getParameter(2);
    fptype sigma = pc.getParameter(1) * (1 + pc.getParameter(3));
    fptype ret   = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    pc.incrementIndex(1, 4, 0, 1, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_ScaledGaussian = device_ScaledGaussian;

__host__ ScaledGaussianPdf::ScaledGaussianPdf(
    std::string n, Observable _x, Variable mean, Variable sigma, Variable delta, Variable epsilon)
    : GooPdf(n, _x) {
    registerParameter(mean);
    registerParameter(sigma);
    registerParameter(delta);
    registerParameter(epsilon);

    initialize();
}

__host__ void ScaledGaussianPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_ScaledGaussian);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_ScaledGaussian");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

} // namespace GooFit
