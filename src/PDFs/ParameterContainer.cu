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

} // namespace GooFit
