#include "goofit/PDFs/combine/ProdPdf.h"
#include "goofit/Log.h"

#include <algorithm>

namespace GooFit {

__device__ fptype device_ProdPdfs(fptype* evt, ParameterContainer &pc) {
    // Index structure is nP | F1 P1 | F2 P2 | ...
    // where nP is number of parameters, Fs are function indices, and Ps are parameter indices

    int numParams = RO_CACHE(pc.parameters[0]);
    fptype ret = 1;

    for(int i = 1; i < numParams; i += 2) {
        pc.incrementIndex (1, numParams + 1, 1, 1, 1);

        //fptype curr = (*(reinterpret_cast<device_function_ptr>(device_function_table[fcnIdx])))(evt, p, paramIndices + parIdx);
        fptype curr = callFunction(evt, pc);

	//TODO: Unsure where this will be located, but needs to be fetched either before or after...
        //curr *= normalisationFactors[parIdx];

        ret *= curr;
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

    observablesList = getObservables(); // Gathers from components

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

__host__ void ProdPdf::recursiveSetIndices () {
    GET_FUNCTION_ADDR(ptr_to_ProdPdfs);
   
    GOOFIT_TRACE("host_function_table[{}] = {}", num_device_functions, getName ());
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;

    populateArrays ();
}

__host__ fptype ProdPdf::normalize() const {
    if(varOverlaps) {
        // Two or more components share an observable and cannot be separately
        // normalized, since \int A*B dx does not equal int A dx * int B dx.
        recursiveSetNormalisation(fptype(1.0));
        //MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);

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

    host_normalisations[normalIdx] = 1;
    //MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);

    return 1.0;
}
} // namespace GooFit
