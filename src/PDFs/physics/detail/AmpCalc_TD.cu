#include <stdio.h>

#include <goofit/PDFs/physics/detail/AmpCalc_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/EvalVar.h>
#include "goofit/utilities/DebugTools.h"
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

struct LS_EntryDebugger_TD 
{
    const unsigned int LS_NUM;
    const Cache_SF_LS_TD_EntryFinder _ENTRY_FINDER;

    LS_EntryDebugger_TD(
        unsigned int lsNum,
        const Cache_SF_LS_TD_EntryFinder& entryFinder)
        : LS_NUM(lsNum),
        _ENTRY_FINDER(entryFinder) {}

    __device__ fpcomplex operator()(unsigned int evtNum) const
    {
        ParameterContainer pc;
        pc.prepare(_ENTRY_FINDER._dalitzFuncId);

        return _ENTRY_FINDER.getLSVal(pc, evtNum, LS_NUM);
    }
};

struct SF_EntryDebugger_TD
{
    const unsigned int SF_NUM;
    const Cache_SF_LS_TD_EntryFinder _ENTRY_FINDER;

    SF_EntryDebugger_TD(
        unsigned int sfNum,
        const Cache_SF_LS_TD_EntryFinder& entryFinder)
        : SF_NUM(sfNum),
        _ENTRY_FINDER(entryFinder) {}

    __device__ fptype operator()(unsigned int evtNum) const
    {
        ParameterContainer pc;
        pc.prepare(_ENTRY_FINDER._dalitzFuncId);

        return _ENTRY_FINDER.getSFVal(pc, evtNum, SF_NUM);
    }
};

AmpCalc_TD::AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx)
    : _entryFinder( Cache_SF_LS_TD_EntryFinder(nPerm, ampIdx) ) {}

__device__ auto AmpCalc_TD::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    //int tid = blockIdx.x*blockDim.x + threadIdx.x;
    //bool PRINT_ME = tid == 0;
    const bool PRINT_ME = true;

    ParameterContainer pc;
    pc.prepare(_entryFinder._dalitzFuncId); 

    const unsigned int LS_STEP = _entryFinder.getLSStep(pc);
    const unsigned int SF_STEP = _entryFinder.getSFStep(pc);

    unsigned int evtNum     = thrust::get<0>(t);

    fpcomplex returnVal(0, 0);

    if(PRINT_ME) 
    {
        printf("_dalitzFuncId: %d\n", _entryFinder._dalitzFuncId);
        printf("_AmpIdx: %i\n", _entryFinder._AMP_IDX);
        printf("_nPerm: %i\n", _entryFinder._N_PERM);
        printf("LS_step: %i\n", LS_STEP);
        printf("SF_step: %i\n", SF_STEP);
        printf("TotAMP: %d\n", _entryFinder.getTotAMP(pc));
        printf("TotLS: %d\n", _entryFinder.getTotLS(pc));
    }

    for(int i = 0; i < _entryFinder._N_PERM; ++i) 
    {
        fpcomplex ret(1, 0);

        if(PRINT_ME) 
        {
            printf("i: %i\n", i);
        }

        for(int j = i * LS_STEP; j < (i + 1) * LS_STEP; ++j)  
        {
            const fpcomplex lsVal = _entryFinder.getLSVal(pc, evtNum, j);
            ret *= lsVal;
            // printf("Lineshape %i = (%.7g, %.7g)\n", j, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP +
            // _AmpIdx + 4 + j]]).real, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 +
            // j]]).imag);
            if(PRINT_ME) 
            {
                printf("LS iteration %i\n", j);
                printf("LS value: (%f, %f)\n", lsVal.real(), lsVal.imag());
                printf("Running ret total: (%f, %f)\n", ret.real(), ret.imag());
            }
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_STEP; j < (i + 1) * SF_STEP; ++j)  
        {
            const fptype sfVal = _entryFinder.getSFVal(pc, evtNum, j);
            ret *= sfVal;
            // printf(" SF = %.7g\n", (cResSF_TD[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 4
            // + numLS + j]].real));
            if(PRINT_ME) 
            {
                printf("SF iteration %i\n", j);
                printf("SF value: %f\n", sfVal);
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
    returnVal *= (1 / sqrt(static_cast<fptype>(_entryFinder._N_PERM)));
    
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real(), returnVal.imag());

    if(PRINT_ME) 
    {
        printf("Value after including nPerm: (%f, %f)\n", returnVal.real(), returnVal.imag());
        printf("Finished AmpCalc (TD).\n\n");
    }
    
    return returnVal;
}

__host__ std::vector<unsigned int> AmpCalc_TD::getLineshapeIndices(int totalAMP) const {
    bool printStatus = false;

    std::vector<unsigned int> hostAmpIndices = DebugTools::copyAmpIndicesToHost();

    unsigned int numLS   = hostAmpIndices[totalAMP + _entryFinder._AMP_IDX];
    unsigned int LS_step = numLS / _entryFinder._N_PERM;

    std::vector<unsigned int> ret;

    if(printStatus) {
        for (const auto& h : hostAmpIndices)
        {
            printf("AmpIdx: %d\n", h);
        }

        printf("In getLineShapesIndices:\n");
        printf("_AmpIdx: %i\n", _entryFinder._AMP_IDX);
        printf("_nPerm: %i\n", _entryFinder._N_PERM);
        printf("totalAMP: %i\n", totalAMP);
        printf("numLS: %i\n", numLS);
        printf("LS_step: %i\n", LS_step);
    }

    for(int i = 0; i < _entryFinder._N_PERM; ++i) {
        if(printStatus) {
            printf("i: %i\n", i);
        }

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            unsigned int idx = hostAmpIndices[totalAMP + _entryFinder._AMP_IDX + 4 + j];
            ret.push_back(idx);

            if(printStatus) {
                printf("LS iteration %i: LS index %i\n", j, idx);
            }
        }
    }

    if(printStatus) 
    {
        for (const auto& l : ret)
        {
            printf("LSIdx: %d\n", l);
        }
        printf("Finished getLineshapeIndices.\n\n");
    }

    return ret;
}

__host__ std::vector<unsigned int> AmpCalc_TD::getSpinFactorIndices(int totalAMP) const {
    bool printStatus = false;

    std::vector<unsigned int> hostAmpIndices = DebugTools::copyAmpIndicesToHost();

    unsigned int numLS   = hostAmpIndices[totalAMP + _entryFinder._AMP_IDX];
    unsigned int numSF   = hostAmpIndices[totalAMP + _entryFinder._AMP_IDX + 1];
    unsigned int SF_step = numSF / _entryFinder._N_PERM;

    std::vector<unsigned int> ret;

    if(printStatus) {
        printf("In getSpinFactorIndices:\n");
        printf("_AmpIdx: %i\n", _entryFinder._AMP_IDX);
        printf("_nPerm: %i\n", _entryFinder._N_PERM);
        printf("totalAMP: %i\n", totalAMP);
        printf("numLS: %i\n", numLS);
        printf("numSF: %i\n", numSF);
        printf("SF_step: %i\n", SF_step);
    }

    for(int i = 0; i < _entryFinder._N_PERM; ++i) {
        if(printStatus) {
            printf("i: %i\n", i);
        }

        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            unsigned int idx = hostAmpIndices[totalAMP + _entryFinder._AMP_IDX + 4 + numLS + j];
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

thrust::host_vector<fpcomplex> AmpCalc_TD::debugLS(
    unsigned int lsNum,
    const thrust::device_vector<unsigned int>& evtNums) const
{
    thrust::device_vector<fpcomplex> results_d(evtNums.size());
    thrust::transform(
        evtNums.cbegin(), 
        evtNums.cend(), 
        results_d.begin(), 
        LS_EntryDebugger_TD(lsNum, _entryFinder));

    thrust::host_vector<fpcomplex> results_h = results_d;
    return results_h;
}

thrust::host_vector<fptype> AmpCalc_TD::debugSF(
    unsigned int sfNum,
    const thrust::device_vector<unsigned int>& evtNums) const
{
    thrust::device_vector<fptype> results_d(evtNums.size());
    thrust::transform(
        evtNums.cbegin(), 
        evtNums.cend(), 
        results_d.begin(), 
        SF_EntryDebugger_TD(sfNum, _entryFinder));

    thrust::host_vector<fptype> results_h = results_d;
    return results_h;
}

} // namespace GooFit
