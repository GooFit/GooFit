#include "goofit/PDFs/combine/ProdPdf.h"
#include "goofit/Log.h"

#include <algorithm>

namespace GooFit {

__device__ fptype device_ProdPdfs(fptype* evt, ParameterContainer &pc) {
    int numCons = RO_CACHE(pc.constants[pc.constantIdx]);
    int numComps = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    fptype ret = 1;

    //pc.incrementIndex (1, 0, numCons, 0, 1);
    pc.incrementIndex();
    for(int i = 0; i < numComps; i ++) {
        fptype norm = pc.normalisations[pc.normalIdx + 1];
        fptype curr = callFunction(evt, pc);

        curr *= norm;

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

        //we push a placeholder that is used to indicate 
        //constantsList.push_back (0);
    }

    observablesList = getObservables(); // Gathers from components

    //Add that we have a components size
    constantsList.push_back(components.size());

    std::vector<Variable *> observableCheck; // Use to check for overlap in observables

    for(PdfBase* p : comps) {
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
   
    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_ProdPdfs");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;

    populateArrays ();
}

__host__ fptype ProdPdf::normalize() const {
    if(varOverlaps) {
        // Two or more components share an observable and cannot be separately
        // normalized, since \int A*B dx does not equal int A dx * int B dx.
        recursiveSetNormalisation(fptype(1.0));
        MEMCPY_TO_SYMBOL(d_normalisations, host_normalisations, totalNormalisations*sizeof(fptype), 0, cudaMemcpyHostToDevice);
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

    host_normalisations[normalIdx + 1] = 1;
    //MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);

    return 1.0;
}
} // namespace GooFit
