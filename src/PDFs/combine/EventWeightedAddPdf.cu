#include <goofit/Error.h>
#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/combine/EventWeightedAddPdf.h>

namespace GooFit {

__device__ auto device_EventWeightedAddPdfs(fptype *evt, ParameterContainer &pc) -> fptype {
    int numConstants = pc.getNumConstants();
    int numObs       = pc.getNumObservables();

    int comps = pc.getConstant(0);

    fptype ret         = 0;
    fptype totalWeight = 0;

    ParameterContainer pci = pc;

    pci.incrementIndex(1, 0, numConstants, numObs, 1);

    for(int i = 0; i < comps - 1; ++i) {
        int id        = pc.getObservable(i);
        fptype norm   = pci.getNormalization(0);
        fptype weight = RO_CACHE(evt[id]);
        totalWeight += weight;
        fptype curr = callFunction(evt, pci);
        ret += weight * curr * norm;
    }

    // numParameters does not count itself. So the array structure for two functions is
    // nP | F P | F P | nO | o1
    // in which nP = 4. and nO = 1. Therefore the parameter index for the last function pointer is nP, and the function
    // index is nP-1.
    // fptype last = (*(reinterpret_cast<device_function_ptr>(d_function_table[indices[numParameters-1]])))(evt, p,
    // paramIndices + indices[numParameters]);

    pc                = pci;
    fptype normFactor = pc.getNormalization(0);

    fptype last = callFunction(evt, pc);
    ret += (1 - totalWeight) * last * normFactor;

    return ret;
}

__device__ auto device_EventWeightedAddPdfsExt(fptype *evt, ParameterContainer &pc) -> fptype {
    // numParameters does not count itself. So the array structure for two functions is
    // nP | F P | F P | nO | o1 o2
    // in which nP = 4, nO = 2.

    int numConstants = pc.getNumConstants();
    int numObs       = pc.getNumObservables();

    int comps = pc.getConstant(0);

    fptype ret         = 0;
    fptype totalWeight = 0;

    ParameterContainer pci = pc;
    pci.incrementIndex(1, 0, numConstants, numObs, 1);

    for(int i = 0; i < comps; ++i) {
        int id        = pc.getObservable(i);
        fptype norm   = pci.getNormalization(0);
        fptype weight = RO_CACHE(evt[id]);
        fptype curr   = callFunction(evt, pci);
        // if ((0 == BLOCKIDX) && (THREADIDX < 5) && (isnan(curr))) printf("NaN component %i %i\n", i, THREADIDX);
        ret += weight * curr * norm;
        totalWeight += weight;

        // if ((gpuDebug & 1) && (0 == THREADIDX))
        // if ((gpuDebug & 1) && (1 > evt[8]))
        // if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
        // printf("EventWeightedExt: %i %f %f | %f %f %f %f %f %f %f\n", i, curr, weight, evt[0], evt[1], evt[2],
        // evt[3], evt[4], evt[5], evt[6]);
        // printf("EventWeightedExt: %i %f %f | %f %f \n", i, curr, weight, normalizationFactors[indices[2*(i+1)]], curr
        // * normalizationFactors[indices[2*(i+1)]]);
        // printf("EventWeightedExt: %i : %i %.10f %.10f %.10f %f %f %f\n", (int) floor(0.5 + evt[8]), i, curr, weight,
        // ret, normalizationFactors[indices[2*(i+1)]], evt[6], evt[7]);
    }

    ret /= totalWeight;

    // update our pc structure.
    pc = pci;

    // if (0 >= ret) printf("Zero sum %f %f %f %f %f %f %f %f %f %f\n", evt[0], evt[1], evt[2], evt[3], evt[4], evt[5],
    // evt[6], evt[7], evt[8], evt[9]);

    return ret;
}

__device__ device_function_ptr ptr_to_EventWeightedAddPdfs    = device_EventWeightedAddPdfs;
__device__ device_function_ptr ptr_to_EventWeightedAddPdfsExt = device_EventWeightedAddPdfsExt;

EventWeightedAddPdf::EventWeightedAddPdf(std::string n, std::vector<Observable> weights, std::vector<PdfBase *> comps)
    : CombinePdf("EventWeightedAddPdf", n) {
    if(weights.size() != comps.size() && (weights.size() + 1) != comps.size())
        throw GooFit::GeneralError("Size of weights {} (+1) != comps {}", weights.size(), comps.size());

    // Indices stores (function index)(function parameter index) doublet for each component.
    // Last component has no weight index unless function is extended. Notice that in this case, unlike
    // AddPdf, weight indices are into the event, not the parameter vector, hence they
    // are not added to the pindices array at this stage, although 'initialize' will reserve space
    // for them.
    for(PdfBase *p : comps) {
        GOOFIT_TRACE("EventWeighted component: {}", p->getName());
        components.push_back(p);
        if(components.back() == nullptr)
            throw GooFit::GeneralError("Invalid component");
    }

    extended = true;

    for(unsigned int w = 0; w < weights.size(); ++w) {
        if(components[w] == nullptr)
            throw GooFit::GeneralError("Invalid component");
        registerObservable(weights[w]);
        // adding room for observable offset
    }

    if(components.back() == nullptr)
        throw GooFit::GeneralError("Invalid component");

    if(weights.size() < components.size()) {
        extended = false;
        // TODO:adding components as parameters to get them to be used (for now)
        // parametersList.push_back(0);
    }

    // This must occur after registering weights, or the indices will be off - the device functions assume that the
    // weights are first.
    observablesList = getObservables();

    registerConstant(components.size());

    if(extended)
        registerFunction("ptr_to_EventWeightedAddPdfsExt", ptr_to_EventWeightedAddPdfsExt);
    else
        registerFunction("ptr_to_EventWeightedAddPdfs", ptr_to_EventWeightedAddPdfs);

    initialize();
}

__host__ auto EventWeightedAddPdf::normalize() -> fptype {
    // if (cpuDebug & 1) std::cout << "Normalizing EventWeightedAddPdf " << getName() << " " << components.size() <<
    // std::endl;

    // Here the PDFs have per-event weights, so there is no per-PDF weight
    // to keep track of. All we can do is normalize the components.
    for(PdfBase *comp : components)
        comp->normalize();

    host_normalizations[normalIdx + 1] = 1.0;
    cachedNormalization                = 1.0;

    return 1.0;
}

} // namespace GooFit
