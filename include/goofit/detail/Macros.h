#pragma once

#include <goofit/Log.h>

// The following macros help with registering new values

#define GOOFIT_START_PDF                                                                                               \
    std::vector<unsigned int> pindices;                                                                                \
    pindices.push_back(0);

#define GOOFIT_ADD_PARAM(i, par, name)                                                                                 \
    {                                                                                                                  \
        pindices.push_back(registerParameter((par)));                                                                  \
        if((i) != pindices.size())                                                                                     \
            throw GooFit::GeneralError(                                                                                \
                "{}: Param {} {} actually has number {}", getName(), (name), (i), pindices.size());                    \
    }

#define GOOFIT_ADD_OBS(par) {registarObservable((par))};

#define GOOFIT_ADD_CONST(i, par, name)                                                                                 \
    {                                                                                                                  \
        pindices.push_back(registerConstants(1));                                                                      \
        MEMCPY_TO_SYMBOL(functorConstants, &(par), sizeof par, cIndex * sizeof par, cudaMemcpyHostToDevice);           \
        if((i) != pindices.size())                                                                                     \
            throw GooFit::GeneralError(                                                                                \
                "{}: Const {} {} actually has number {}", getName(), (name), (i), pindices.size());                    \
        GOOFIT_DEBUG("{}: Registered constant value {}={} at c:{}", getName(), (name), (par), pindices.size());        \
    }

#define GOOFIT_ADD_INT(i, par, name)                                                                                   \
    {                                                                                                                  \
        GOOFIT_DEBUG("{}: Registered integer {} at {}", getName(), (name), (par), pindices.size());                    \
        pindices.push_back((par));                                                                                     \
        if((i) != pindices.size())                                                                                     \
            throw GooFit::GeneralError(                                                                                \
                "{}: Int {} {} actually has number {}", getName(), (name), (i), pindices.size());                      \
    }

#define GOOFIT_FINALIZE_PDF                                                                                            \
    initialize(pindices);                                                                                              \
    GOOFIT_DEBUG("{}: Initializing indices", getName());

#define GOOFIT_GET_PARAM(i) cudaArray[indices[(i)]]
#define GOOFIT_GET_INT(i) indices[(i)]
#define GOOFIT_GET_CONST(i) functorConstants[indices[(i)]]

#define GOOFIT_PDF_IMPL_1(n) __device__ fptype #n(fptype *evt, fptype *p, unsigned int *indices)
#define GOOFIT_PDF_IMPL_3(n) __device__ fpcomplex #n(fptype Mpair, fptype m1, fptype m2, unsigned int *indices)
