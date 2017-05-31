#include "goofit/PDFs/combine/ProdPdf.h"
#include <algorithm>

namespace GooFit {

__device__ fptype device_ProdPdfs(fptype *evt, fptype *p, unsigned int *indices) {
    // Index structure is nP | F1 P1 | F2 P2 | ...
    // where nP is number of parameters, Fs are function indices, and Ps are parameter indices

    int numParams = RO_CACHE(indices[0]);
    fptype ret    = 1;

    for(int i = 1; i < numParams; i += 2) {
        int fcnIdx = RO_CACHE(indices[i + 0]);
        int parIdx = RO_CACHE(indices[i + 1]);

        // fptype curr = (*(reinterpret_cast<device_function_ptr>(device_function_table[fcnIdx])))(evt, p, paramIndices
        // + parIdx);
        fptype curr = callFunction(evt, fcnIdx, parIdx);
        curr *= normalisationFactors[parIdx];
        // if ((isnan(ret)) || (isnan(curr)) || (isnan(normalisationFactors[parIdx])) || (isinf(ret)) || (isinf(curr)))
        // printf("device_Prod 2: (%f %f %f %f %f) %f %f %f %i %i %i\n", evt[0], evt[1], evt[2], evt[3], evt[4], curr,
        // ret, normalisationFactors[parIdx], i, parIdx, numParams);
        ret *= curr;

        // if ((0 == THREADIDX) && (0 == BLOCKIDX) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // if (0.0001 < ret)
        // if ((gpuDebug & 1) && (isnan(curr)) && (paramIndices + debugParamIndex == indices))
        // if ((isnan(ret)) || (isnan(curr)) || (isnan(normalisationFactors[parIdx])))
        // printf("device_Prod: (%f %f %f %f %f) %f %f %f %i %i %i\n", evt[0], evt[1], evt[2], evt[3], evt[4], curr,
        // ret, normalisationFactors[parIdx], i, parIdx, numParams);
        // printf("(%i, %i) device_Prod: (%f %f %f %f) %f %f %f %i\n", BLOCKIDX, THREADIDX, evt[0], evt[8], evt[6],
        // evt[7], curr, ret, normalisationFactors[parIdx], i);
        // printf("(%i, %i) device_Prod: (%f %f) %f %f %f %i\n", BLOCKIDX, THREADIDX, evt[0], evt[1], curr, ret,
        // normalisationFactors[parIdx], i);
    }

    return ret;
}

__device__ device_function_ptr ptr_to_ProdPdfs = device_ProdPdfs;

ProdPdf::ProdPdf(std::string n, std::vector<PdfBase *> comps)
    : GooPdf(nullptr, n)
    , varOverlaps(false) {
    std::vector<unsigned int> pindices;

    for(PdfBase *p : comps) {
        components.push_back(p);
    }

    observables = getObservables(); // Gathers from components

    std::vector<Variable *> observableCheck; // Use to check for overlap in observables

    // Indices stores (function index)(function parameter index)(variable index) for each component.
    for(PdfBase *p : comps) {
        pindices.push_back(p->getFunctionIndex());
        pindices.push_back(p->getParameterIndex());

        if(varOverlaps)
            continue; // Only need to establish this once.

        std::vector<Variable *> currObses = p->getObservables();

        for(Variable *o : currObses) {
            if(find(observableCheck.begin(), observableCheck.end(), o) == observableCheck.end())
                continue;

            varOverlaps = true;
            break;
        }

        observableCheck = p->getObservables();
    }

    if(varOverlaps) { // Check for components forcing separate normalisation
        for(PdfBase *p : comps) {
            if(p->getSpecialMask() & PdfBase::ForceSeparateNorm)
                varOverlaps = false;
        }
    }

    GET_FUNCTION_ADDR(ptr_to_ProdPdfs);
    initialize(pindices);
}

__host__ fptype ProdPdf::normalize() const {
    if(varOverlaps) {
        // Two or more components share an observable and cannot be separately
        // normalized, since \int A*B dx does not equal int A dx * int B dx.
        recursiveSetNormalisation(fptype(1.0));
        MEMCPY_TO_SYMBOL(
            normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);

        // Normalize numerically.
        // std::cout << "Numerical normalisation of " << getName() << " due to varOverlaps.\n";
        fptype ret = GooPdf::normalize();
        // if (cpuDebug & 1)
        // std::cout << "ProdPdf " << getName() << " has normalisation " << ret << " " << host_callnumber << std::endl;
        return ret;
    }

    // Normalize components individually
    for(PdfBase *c : components) {
        c->normalize();
    }

    host_normalisation[parameters] = 1;
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    return 1.0;
}
} // namespace GooFit
