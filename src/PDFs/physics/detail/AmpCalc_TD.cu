#include <goofit/PDFs/physics/detail/AmpCalc_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/EvalVar.h>
#include "goofit/utilities/DebugTools.h"
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

AmpCalc_TD::AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx, unsigned int CacheIdx)
    : _nPerm(nPerm)
    , _AmpIdx(ampIdx)
    , _CacheIdx(CacheIdx) {}

__device__ auto AmpCalc_TD::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int cacheToUse = pc.getConstant(6);
    unsigned int totalLS    = pc.getConstant(7);
    unsigned int totalSF    = pc.getConstant(8);
    unsigned int totalAMP   = pc.getConstant(9);
    unsigned int offset     = totalLS + totalSF;
    unsigned int numLS      = AmpIndices[_CacheIdx][totalAMP + _AmpIdx];
    unsigned int numSF      = AmpIndices[_CacheIdx][totalAMP + _AmpIdx + 1];
    unsigned int evtNum     = thrust::get<0>(t);

    fpcomplex returnVal(0, 0);
    unsigned int SF_step = numSF / _nPerm;
    unsigned int LS_step = numLS / _nPerm;

    // int tid = blockIdx.x*blockDim.x + threadIdx.x;
    // bool printStatus = tid == 0;
    bool printStatus = false;

    if(printStatus) {
        printf("_AmpIdx: %i\n", _AmpIdx);
        printf("_nPerm: %i\n", _nPerm);
        printf("totalAMP: %i\n", totalAMP);
        printf("numLS: %i\n", numLS);
        printf("LS_step: %i\n", LS_step);
        printf("numSF: %i\n", numSF);
        printf("SF_step: %i\n", SF_step);
    }

    for(int i = 0; i < _nPerm; ++i) {
        fpcomplex ret(1, 0);

        if(printStatus) {
            printf("i: %i\n", i);
        }

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            int idx = AmpIndices[cacheToUse][totalAMP + _AmpIdx + 4 + j];
            ret *= (cResSF_TD[cacheToUse][evtNum * offset + idx]);
            // printf("Lineshape %i = (%.7g, %.7g)\n", j, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP +
            // _AmpIdx + 4 + j]]).real, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 +
            // j]]).imag);
            if(printStatus) {
                printf("LS iteration %i: LS index %i\n", j, idx);
            }
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            int idx = AmpIndices[cacheToUse][totalAMP + _AmpIdx + 4 + numLS + j];
            ret *= (cResSF_TD[cacheToUse][evtNum * offset + totalLS + idx].real());
            // printf(" SF = %.7g\n", (cResSF_TD[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 4
            // + numLS + j]].real));
            if(printStatus) {
                printf("SF iteration %i: SF index %i\n", j, idx);
            }
        }

        // printf("Lineshape Product * SF = (%.7g, %.7g)\n", ret.real, ret.imag);

        returnVal += ret;
    }

    if(printStatus) {
        printf("Finished AmpCalc (TD).\n\n");
    }

    returnVal *= (1 / sqrt(static_cast<fptype>(_nPerm)));
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real(), returnVal.imag());
    return returnVal;
}

__host__ std::vector<unsigned int> AmpCalc_TD::getLineshapeIndices(int totalAMP) const {
    bool printStatus = false;

    std::vector<unsigned int> hostAmpIndices = DebugTools::copyAmpIndicesToHost();
    //printf("_CacheIdx: %i\n");
    //printf("hostAmpIndices\n");
    //for(auto idx:hostAmpIndices){
    //    printf("%i ",idx);
    //}
    //printf("\n");
    unsigned int numLS   = hostAmpIndices[(_CacheIdx * 500) + totalAMP + _AmpIdx]; //offset the indices by the correct number of Amp4Body_TD PDFs there
    unsigned int LS_step = numLS / _nPerm;

    std::vector<unsigned int> ret;

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
            unsigned int idx = hostAmpIndices[_CacheIdx * 500 + totalAMP + _AmpIdx + 4 + j];
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

__host__ std::vector<unsigned int> AmpCalc_TD::getSpinFactorIndices(int totalAMP) const {
    bool printStatus = false;

    std::vector<unsigned int> hostAmpIndices = DebugTools::copyAmpIndicesToHost();

    unsigned int numLS   = hostAmpIndices[_CacheIdx * 500 + totalAMP + _AmpIdx];
    unsigned int numSF   = hostAmpIndices[_CacheIdx * 500 + totalAMP + _AmpIdx + 1];
    unsigned int SF_step = numSF / _nPerm;

    std::vector<unsigned int> ret;

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
            unsigned int idx = hostAmpIndices[_CacheIdx * 500 + totalAMP + _AmpIdx + 4 + numLS + j];
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
