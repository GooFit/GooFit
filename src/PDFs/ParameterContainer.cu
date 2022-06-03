#include <goofit/PDFs/ParameterContainer.h>

namespace GooFit {

__device__ ParameterContainer::ParameterContainer()
    : parameters(d_parameters)
    , constants(d_constants)
    , observables(d_observables)
    , normalizations(d_normalizations) {}

__device__ ParameterContainer::ParameterContainer(const ParameterContainer &pc)
    : ParameterContainer() {
    parameterIdx  = pc.parameterIdx;
    constantIdx   = pc.constantIdx;
    observableIdx = pc.observableIdx;
    normalIdx     = pc.normalIdx;
    funcIdx       = pc.funcIdx;
}

/// Only on CPU and with Trace enabled, this will run slowly but check to make sure you are correctly incrementing
/// indices
__device__ void
ParameterContainer::incrementIndex(const int funcs, const int params, const int cons, const int obs, const int norms) {
#if !defined(__CUDACC__) && defined(GOOFIT_TRACE_FLAG)
    if(funcs != 1)
        throw GeneralError("Haven't got a clue on how to proceed with incrementIndex checking, sorry");
    if(parameters[parameterIdx] != params || constants[constantIdx] != cons || observables[observableIdx] != obs
       || normalizations[normalIdx] != norms)
        throw GeneralError(
            "Wrong parameters given to incrementIndex(1, {}, {}, {}, {}), should have been (1, {}, {}, {}, {})",
            params,
            cons,
            obs,
            norms,
            parameters[parameterIdx],
            constants[constantIdx],
            observables[observableIdx],
            normalizations[normalIdx]);
#endif

    funcIdx += funcs;
    parameterIdx += params + 1;
    constantIdx += cons + 1;
    observableIdx += obs + 1;
    normalIdx += norms + 1;
}

__device__ void ParameterContainer::incrementIndex() {
    funcIdx++;

    int np = RO_CACHE(parameters[parameterIdx]);
    int nc = RO_CACHE(constants[constantIdx]);
    int no = RO_CACHE(observables[observableIdx]);
    int nn = RO_CACHE(normalizations[normalIdx]);

    parameterIdx += np + 1;
    constantIdx += nc + 1;
    observableIdx += no + 1;
    normalIdx += nn + 1;
}

} // namespace GooFit
