#include <stdio.h>

#include <goofit/PDFs/physics/detail/AmpCalc_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/EvalVar.h>
#include "goofit/utilities/DebugTools.h"
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

AmpCalc_TD::AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx)
    : _nPerm(nPerm)
    , _AmpIdx(ampIdx) {}

__device__ auto AmpCalc_TD::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    // int tid = blockIdx.x*blockDim.x + threadIdx.x;
    // bool printStatus = tid == 0;
    const bool PRINT_ME = false;

    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int cacheToUse = pc.getConstant(5);
    unsigned int totalLS    = pc.getConstant(6);
    unsigned int totalSF    = pc.getConstant(7);
    unsigned int totalAMP   = pc.getConstant(8);
    unsigned int offset     = totalLS + totalSF;
    unsigned int numLS      = AmpIndices[totalAMP + _AmpIdx];
    unsigned int numSF      = AmpIndices[totalAMP + _AmpIdx + 1];
    unsigned int evtNum     = thrust::get<0>(t);

    fpcomplex returnVal(0, 0);
    unsigned int SF_step = numSF / _nPerm;
    unsigned int LS_step = numLS / _nPerm;

    if(PRINT_ME) 
    {
        printf("_AmpIdx: %i\n", _AmpIdx);
        printf("_nPerm: %i\n", _nPerm);
        printf("totalAMP: %i\n", totalAMP);
        printf("numLS: %i\n", numLS);
        printf("LS_step: %i\n", LS_step);
        printf("numSF: %i\n", numSF);
        printf("SF_step: %i\n", SF_step);
    }

    for(int i = 0; i < _nPerm; ++i) 
    {
        fpcomplex ret(1, 0);

        if(PRINT_ME) 
        {
            printf("i: %i\n", i);
        }

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) 
        {
            int idx = AmpIndices[totalAMP + _AmpIdx + 4 + j];
            fpcomplex lsVal = (cResSF_TD[cacheToUse][evtNum * offset + idx]);
            ret *= lsVal;
            // printf("Lineshape %i = (%.7g, %.7g)\n", j, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP +
            // _AmpIdx + 4 + j]]).real, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 +
            // j]]).imag);
            if(PRINT_ME) 
            {
                printf("LS iteration %i: LS index %i\n", j, idx);
                printf("LS value: (%f, %f)\n", lsVal.real(), lsVal.imag());
                printf("Running ret total: (%f, %f)\n", ret.real(), ret.imag());
            }
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) 
        {
            int idx = AmpIndices[totalAMP + _AmpIdx + 4 + numLS + j];
            fpcomplex sfVal = (cResSF_TD[cacheToUse][evtNum * offset + totalLS + idx].real());
            ret *= sfVal;
            // printf(" SF = %.7g\n", (cResSF_TD[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 4
            // + numLS + j]].real));
            if(PRINT_ME) 
            {
                printf("SF iteration %i: SF index %i\n", j, idx);
                printf("SF value: (%f, %f)\n", sfVal.real(), sfVal.imag());
                printf("Running ret total: (%f, %f)\n", ret.real(), ret.imag());
            }
        }

        // printf("Lineshape Product * SF = (%.7g, %.7g)\n", ret.real, ret.imag);

        returnVal += ret;
        if (PRINT_ME)
        {
            printf("Running returnVal total: (%f, %f)\n", returnVal.real(), returnVal.imag());
        }
    }

    if (PRINT_ME)
    {
        printf("Value before including nPerm: (%f, %f)\n", returnVal.real(), returnVal.imag());
    }
    returnVal *= (1 / sqrt(static_cast<fptype>(_nPerm)));
    
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real(), returnVal.imag());

    if(PRINT_ME) 
    {
        printf("Value after including nPerm: (%f, %f)\n", returnVal.real(), returnVal.imag());
        printf("Finished AmpCalc (TD).\n\n");
    }
    
    return returnVal;
}

__host__ std::vector<size_t> AmpCalc_TD::getLineshapeIndices(int totalAMP) const {
    bool printStatus = false;

    std::vector<size_t> hostAmpIndices = DebugTools::copyAmpIndicesToHost();

    unsigned int numLS   = hostAmpIndices[totalAMP + _AmpIdx];
    unsigned int LS_step = numLS / _nPerm;

    std::vector<size_t> ret;

    if(printStatus) {
        printf("In getLineShapesIndices:\n");
        printf("_AmpIdx: %i\n", _AmpIdx);
        printf("_nPerm: %i\n", _nPerm);
        printf("totalAMP: %i\n", totalAMP);
        printf("numLS: %i\n", numLS);
        printf("LS_step: %i\n", LS_step);
    }

    for(int i = 0; i < _nPerm; ++i) {
        if(printStatus) {
            printf("i: %i\n", i);
        }

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            unsigned int idx = hostAmpIndices[totalAMP + _AmpIdx + 4 + j];
            ret.push_back(idx);

            if(printStatus) {
                printf("LS iteration %i: LS index %i\n", j, idx);
            }
        }
    }

    if(printStatus) {
        printf("Finished getLineshapeIndices.\n\n");
    }

    return ret;
}

__host__ std::vector<size_t> AmpCalc_TD::getSpinFactorIndices(int totalAMP) const {
    bool printStatus = false;

    std::vector<size_t> hostAmpIndices = DebugTools::copyAmpIndicesToHost();

    unsigned int numLS   = hostAmpIndices[totalAMP + _AmpIdx];
    unsigned int numSF   = hostAmpIndices[totalAMP + _AmpIdx + 1];
    unsigned int SF_step = numSF / _nPerm;

    std::vector<size_t> ret;

    if(printStatus) {
        printf("In getSpinFactorIndices:\n");
        printf("_AmpIdx: %i\n", _AmpIdx);
        printf("_nPerm: %i\n", _nPerm);
        printf("totalAMP: %i\n", totalAMP);
        printf("numLS: %i\n", numLS);
        printf("numSF: %i\n", numSF);
        printf("SF_step: %i\n", SF_step);
    }

    for(int i = 0; i < _nPerm; ++i) {
        if(printStatus) {
            printf("i: %i\n", i);
        }

        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            unsigned int idx = hostAmpIndices[totalAMP + _AmpIdx + 4 + numLS + j];
            ret.push_back(idx);

            if(printStatus) {
                printf("SF iteration %i: SF index %i\n", j, idx);
            }
        }
    }

    if(printStatus) {
        printf("Finished getSpinFactorIndices.\n\n");
    }

    return ret;
}

} // namespace GooFit
