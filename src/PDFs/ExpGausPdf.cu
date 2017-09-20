#include "goofit/PDFs/basic/ExpGausPdf.h"

namespace GooFit {

__device__ fptype device_ExpGaus(fptype *evt, ParameterContainer &pc) {
    int id = pc.constants[pc.constantIdx + 2];

    fptype x     = evt[id];
    fptype mean  = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype sigma = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype alpha = RO_CACHE(pc.parameters[pc.parameterIdx + 3]);

    fptype ret    = 0.5 * alpha;
    fptype exparg = ret * (2 * mean + alpha * sigma * sigma - 2 * x);
    fptype erfarg = (mean + alpha * sigma * sigma - x) / (sigma * 1.4142135623);

    ret *= exp(exparg);
    ret *= erfc(erfarg);

    pc.incrementIndex (1, 3, 2, 0, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_ExpGaus = device_ExpGaus;

ExpGausPdf::ExpGausPdf(std::string n, Variable *_x, Variable *mean, Variable *sigma, Variable *tau)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;

    //reserving index for _x    
    constantsList.push_back (observablesList.size());
    constantsList.push_back (0);

    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigma));
    pindices.push_back(registerParameter(tau));
    GET_FUNCTION_ADDR(ptr_to_ExpGaus);
    initialize(pindices);
}

__host__ void ExpGausPdf::recursiveSetIndices () {
    GET_FUNCTION_ADDR(ptr_to_ExpGaus);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_ExpGaus");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;

    populateArrays ();
 
}

} // namespace GooFit
