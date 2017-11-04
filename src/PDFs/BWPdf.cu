#include "goofit/PDFs/basic/BWPdf.h"

namespace GooFit {

__device__ fptype device_BW(fptype *evt, ParameterContainer &pc) {
    int id = RO_CACHE(pc.observables[pc.observableIdx + 1]);

    fptype x      = evt[id];
    fptype mean   = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype gamma  = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype rootPi = -2. * atan2(-1.0, 0.0);
    fptype ret    = (gamma / ((x - mean) * (x - mean) + gamma * gamma / 4)) / (2 * rootPi);

    pc.incrementIndex(1, 2, 0, 1, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_BW = device_BW;

__host__ BWPdf::BWPdf(std::string n, Variable *_x, Variable *mean, Variable *width)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;

    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(width));
    GET_FUNCTION_ADDR(ptr_to_BW);
    initialize(pindices);
}

__host__ void BWPdf::recursiveSetIndices() {
    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_BW");
    GET_FUNCTION_ADDR(ptr_to_BW);

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

} // namespace GooFit
