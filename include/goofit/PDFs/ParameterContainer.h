#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/// This is a container that is used to communicate to the device PDF functions
struct ParameterContainer {
    __host__ __device__ ParameterContainer();
    __host__ __device__ ParameterContainer(const ParameterContainer &pc);

    fptype *parameters;
    fptype *constants;
    fptype *observables;
    fptype *normalizations;

    int parameterIdx{0};
    int constantIdx{0};
    int observableIdx{0};
    int normalIdx{0};

    int funcIdx{0};

    inline __device__ fptype getParameter(const int i) { return RO_CACHE(parameters[parameterIdx + i + 1]); }

    inline __device__ fptype getConstant(const int i) { return RO_CACHE(constants[constantIdx + i + 1]); }

    inline __device__ fptype getObservable(const int i) { return RO_CACHE(observables[observableIdx + i + 1]); }

    inline __device__ fptype getNormalization(const int i) { return RO_CACHE(normalizations[normalIdx + i + 1]); }

    inline __device__ int getNumParameters() { return (int)RO_CACHE(parameters[parameterIdx]); }

    inline __device__ int getNumConstants() { return (int)RO_CACHE(constants[constantIdx]); }

    inline __device__ int getNumObservables() { return (int)RO_CACHE(observables[observableIdx]); }

    inline __device__ int getNumNormalizations() { return (int)RO_CACHE(normalizations[normalIdx]); }

    // each PDF needs to supply the amount of each array used.
    // This function automatically adds +1 for the size.
    __device__ void incrementIndex(const int funcs, const int params, const int cons, const int obs, const int norms);

    // slow version, avoid at all costs!
    __device__ void incrementIndex();
};

} // namespace GooFit
