#include "goofit/PDFs/basic/GaussianPdf.h"
#include "goofit/Log.h"

namespace GooFit {



__device__ fptype device_Gaussian(fptype* evt, ParameterContainer &pc) {
    fptype x = evt[0];
    fptype mean = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype sigma = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);

    pc.incrementIndex (1, 2, 0, 0, 1);
    fptype ret = exp(-0.5*(x-mean)*(x-mean)/(sigma*sigma));

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

__host__ void GaussianPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_Gaussian);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_Gaussian");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;

    populateArrays ();
}

__host__ fptype GaussianPdf::integrate(fptype lo, fptype hi) const {
    static const fptype rootPi = sqrt(atan2(0.0, -1.0));

    // Integral over all R.
    fptype sigma = host_parameters[parametersIdx + 2];
    sigma *= root2*rootPi;
    return sigma;
}

} // namespace GooFit
