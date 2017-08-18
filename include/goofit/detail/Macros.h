#pragma once

// The following macros help with registering new values


#define GOOFIT_ADD_PARAM(pindices,par) \
    (pindices).push_back(registerParameter((par));

#define GOOFIT_ADD_OBS(par) \
    registarObservable((par));

#define GOOFIT_ADD_CONST(pindices, par, name) \
            (pindices).push_back(registerConstants); \
            MEMCPY_TO_SYMBOL(functorConstants, &(par), sizeof par, cIndex*sizeof par, cudaMemcpyHostToDevice); \
            GOOFIT_DEBUG("{}: Registered constant value {}={} at c:{}", getName(), (name), (par), cIndex); \

#define GOOFIT_ADD_INT(pindices, par, name) \
            GOOFIT_DEBUG("{}: Registered integer {} at {}", getName(), (name),  (par), (pindices).size()); \
            (pindices).push_back((par));

#define GOOFIT_FINALIZE_PARAM(pindices) \
    initialize((pindices)); \
    GOOFIT_DEBUG("{}: Initializing indices", getName());
