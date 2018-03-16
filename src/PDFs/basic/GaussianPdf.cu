#include <goofit/Log.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/ParameterContainer.h>

namespace GooFit {

__device__ fptype device_Gaussian(fptype *evt, ParameterContainer &pc) {
    int id       = pc.getObservable(0);
    fptype mean  = pc.getParameter(0);
    fptype sigma = pc.getParameter(1);
    fptype x     = evt[id];

    pc.incrementIndex(1, 2, 0, 1, 1);
    fptype ret = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    return ret;
}

__device__ device_function_ptr ptr_to_Gaussian = device_Gaussian;

__host__ GaussianPdf::GaussianPdf(std::string n, Observable _x, Variable mean, Variable sigma)
    : GooPdf(n, _x) {
    registerParameter(mean);
    registerParameter(sigma);

    initialize();
}

__host__ void GaussianPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_Gaussian);
    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_Gaussian");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

__host__ fptype GaussianPdf::integrate(fptype lo, fptype hi) const {
    static const fptype rootPi = sqrt(atan2(0.0, -1.0));

    // Integral over all R.
    fptype sigma = host_parameters[parametersIdx + 2];
    sigma *= root2 * rootPi;
    return sigma;
}

} // namespace GooFit
