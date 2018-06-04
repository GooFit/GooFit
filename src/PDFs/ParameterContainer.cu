#include <goofit/PDFs/ParameterContainer.h>

namespace GooFit {

__device__ __host__ ParameterContainer::ParameterContainer()
    : parameters(d_parameters)
    , constants(d_constants)
    , observables(d_observables)
    , normalisations(d_normalisations) {}

__device__ __host__ ParameterContainer::ParameterContainer(const ParameterContainer &pc)
    : parameters(d_parameters)
    , constants(d_constants)
    , observables(d_observables)
    , normalisations(d_normalisations) {
    parameterIdx  = pc.parameterIdx;
    constantIdx   = pc.constantIdx;
    observableIdx = pc.observableIdx;
    normalIdx     = pc.normalIdx;
    funcIdx       = pc.funcIdx;
}

__device__ void
ParameterContainer::incrementIndex(const int funcs, const int params, const int cons, const int obs, const int norms) {
    // Only on CPU, slow
#if !defined(__CUDACC__) && defined(GOOFIT_TRACE_FLAG)
    if(funcs != 1)
        throw GeneralError("Haven't got a clue on how to procede with incrementIndex checking, sorry");
    if(parameters[parameterIdx] != params || constants[constantIdx] != cons || observables[observableIdx] != obs
       || normalisations[normalIdx] != norms)
        throw GeneralError(
            "Wrong parameters given to incrementIndex(1, {}, {}, {}, {}), should have been (1, {}, {}, {}, {})",
            params,
            cons,
            obs,
            norms,
            parameters[parameterIdx],
            constants[constantIdx],
            observables[observableIdx],
            normalisations[normalIdx]);
#endif

    funcIdx += funcs;
    parameterIdx += params + 1;
    constantIdx += cons + 1;
    observableIdx += obs + 1;
    normalIdx += norms + 1;
}

__device__ void ParameterContainer::incrementIndex() {
    funcIdx++;

    int np = parameters[parameterIdx];
    int nc = constants[constantIdx];
    int no = observables[observableIdx];
    int nn = normalisations[normalIdx];

    parameterIdx += np + 1;
    constantIdx += nc + 1;
    observableIdx += no + 1;
    normalIdx += nn + 1;
}

} // namespace GooFit
