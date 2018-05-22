#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/StepPdf.h>

namespace GooFit {

__device__ fptype device_Step(fptype *evt, ParameterContainer &pc) {
    int id    = pc.getObservable(0);
    fptype x  = evt[id];
    fptype x0 = pc.getParameter(0);
    pc.incrementIndex(1, 1, 0, 1, 1);
    return (x > x0 ? 1 : 0);
}

__device__ device_function_ptr ptr_to_Step = device_Step;
device_function_ptr hptr_to_Step           = device_Step;

__host__ StepPdf::StepPdf(std::string n, Observable _x, Variable x0)
    : GooPdf(n, _x) {
    registerParameter(x0);

    initialize();
}

__host__ fptype StepPdf::integrate(fptype lo, fptype hi) const {
    // unsigned int *indices = host_indices + parameters;
    fptype x0 = parametersList[0].getValue();
    return (hi - x0);
}

__host__ void StepPdf::recursiveSetIndices() { GOOFIT_RECURSIVE_SET_INDICIES(ptr_to_Step); }

} // namespace GooFit
