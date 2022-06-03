#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/// This is a container that is used to communicate to the device PDF functions
struct ParameterContainer {
    __device__ ParameterContainer();
    __device__ ParameterContainer(const ParameterContainer &pc);

    fptype *parameters;
    fptype *constants;
    fptype *observables;
    fptype *normalizations;

    int parameterIdx{0};
    int constantIdx{0};
    int observableIdx{0};
    int normalIdx{0};

    int funcIdx{0};

    inline __device__ fptype getParameter(const int i) const { return RO_CACHE(parameters[parameterIdx + i + 1]); }

    inline __device__ fptype getConstant(const int i) const { return RO_CACHE(constants[constantIdx + i + 1]); }

    inline __device__ fptype getObservable(const int i) const { return RO_CACHE(observables[observableIdx + i + 1]); }

    inline __device__ fptype getNormalization(const int i) const { return RO_CACHE(normalizations[normalIdx + i + 1]); }

    inline __device__ int getNumParameters() const { return (int)RO_CACHE(parameters[parameterIdx]); }

    inline __device__ int getNumConstants() const { return (int)RO_CACHE(constants[constantIdx]); }

    inline __device__ int getNumObservables() const { return (int)RO_CACHE(observables[observableIdx]); }

    inline __device__ int getNumNormalizations() const { return (int)RO_CACHE(normalizations[normalIdx]); }

    // each PDF needs to supply the amount of each array used.
    // This function automatically adds +1 for the size.
    // __device__ void incrementIndex(const int funcs, const int params, const int cons, const int obs, const int norms);

    inline __device__ void
    incrementIndex(const int funcs, const int params, const int cons, const int obs, const int norms) {
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


    // slow version, avoid at all costs!
    inline __device__ void incrementIndex() {
        funcIdx++;

        int np = parameters[parameterIdx];
        int nc = constants[constantIdx];
        int no = observables[observableIdx];
        int nn = normalizations[normalIdx];

        parameterIdx += np + 1;
        constantIdx += nc + 1;
        observableIdx += no + 1;
        normalIdx += nn + 1;
    }
};

} // namespace GooFit
