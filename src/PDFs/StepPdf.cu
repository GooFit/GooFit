#include "goofit/PDFs/basic/StepPdf.h"

namespace GooFit {

__device__ fptype device_Step(fptype *evt, ParameterContainer &pc) {
    int id = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    fptype x  = evt[id];
    fptype x0 = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    pc.incrementIndex (1, 1, 0, 1, 1);
    return (x > x0 ? 1 : 0);
}

__device__ device_function_ptr ptr_to_Step = device_Step;
device_function_ptr hptr_to_Step           = device_Step;

__host__ StepPdf::StepPdf(std::string n, Variable *_x, Variable *x0)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(x0));
    GET_FUNCTION_ADDR(ptr_to_Step);
    initialize(pindices);
}

__host__ fptype StepPdf::integrate(fptype lo, fptype hi) const {
    //unsigned int *indices = host_indices + parameters;
    fptype x0             = parametersList[0]->getValue();
    return (hi - x0);
}

__host__ void StepPdf::recursiveSetIndices () {
    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_Step");
    GET_FUNCTION_ADDR(ptr_to_Step);

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions ++;

    populateArrays ();
}

} // namespace GooFit
