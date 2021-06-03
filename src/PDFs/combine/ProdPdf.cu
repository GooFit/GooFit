#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/combine/ProdPdf.h>

#include <algorithm>

namespace GooFit {

__device__ auto device_ProdPdfs(fptype *evt, ParameterContainer &pc) -> fptype {
    int numCons  = pc.getNumConstants();
    int numComps = pc.getConstant(0);
    int numObs   = pc.getNumObservables();
    fptype ret   = 1;

    pc.incrementIndex(1, 0, numCons, numObs, 1);
    // pc.incrementIndex();
    for(int i = 0; i < numComps; i++) {
        fptype norm = pc.getNormalization(0);
        fptype curr = callFunction(evt, pc);

        curr *= norm;
        ret *= curr;
    }

    return ret;
}

__device__ device_function_ptr ptr_to_ProdPdfs = device_ProdPdfs;

ProdPdf::ProdPdf(std::string n, std::vector<PdfBase *> comps)
    : CombinePdf("ProdPdf", n)
    , varOverlaps(false) {
    for(PdfBase *p : comps) {
        components.push_back(p);
        // we push a placeholder that is used to indicate
        // constantsList.push_back (0);
    }

    observablesList = getObservables(); // Gathers from components

    // Add that we have a components size
    registerConstant(components.size());

    std::vector<Observable> observableCheck; // Use to check for overlap in observables

    for(PdfBase *p : comps) {
        if(varOverlaps)
            continue; // Only need to establish this once.

        std::vector<Observable> currObses = p->getObservables();

        for(Observable &o : currObses) {
            if(find(observableCheck.begin(), observableCheck.end(), o) == observableCheck.end())
                continue;

            varOverlaps = true;
            break;
        }

        observableCheck = p->getObservables();
    }

    if(varOverlaps) { // Check for components forcing separate normalization
        for(PdfBase *p : comps) {
            if(p->getSeparateNorm())
                varOverlaps = false;
        }
    }

    registerFunction("ptr_to_ProdPdfs", ptr_to_ProdPdfs);

    initialize();
}

__host__ auto ProdPdf::normalize() -> fptype {
    if(varOverlaps) {
        // Two or more components share an observable and cannot be separately
        // normalized, since \int A*B dx does not equal int A dx * int B dx.
        recursiveSetNormalization(1.0);
        host_normalizations.sync(d_normalizations);

        // Normalize numerically.
        // std::cout << "Numerical normalization of " << getName() << " due to varOverlaps.\n";
        fptype ret = GooPdf::normalize();
        // if (cpuDebug & 1)
        // std::cout << "ProdPdf " << getName() << " has normalization " << ret << " " << host_callnumber << std::endl;
        return ret;
    }

    // Normalize components individually
    for(PdfBase *c : components) {
        c->normalize();
    }

    host_normalizations.at(normalIdx + 1) = 1.0;
    cachedNormalization                   = 1.0;

    return 1.0;
}
} // namespace GooFit
