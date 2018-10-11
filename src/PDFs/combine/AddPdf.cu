#include <goofit/Error.h>
#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/detail/ThrustOverride.h>

#include <thrust/iterator/constant_iterator.h>
#include <thrust/transform_reduce.h>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

__device__ fptype device_AddPdfs(fptype *evt, ParameterContainer &pc) {
    int numParameters  = pc.getNumParameters();
    fptype ret         = 0;
    fptype totalWeight = 0;

    // make a copy of our parameter container
    ParameterContainer pci = pc;

    // We only call increment once we read our weight/norm for the first iteration.
    pci.incrementIndex();

    for(int i = 0; i < numParameters; i++) {
        // fetch our values from AddPdf
        fptype weight = pc.getParameter(i);
        totalWeight += weight;

        // This is the normal value for the 'callFunction' PDF, so we read from pci
        fptype norm = pci.getNormalization(0);

        // call the first function to add in our PDF.
        fptype curr = callFunction(evt, pci);

        ret += weight * curr * norm;
    }

    // restore our new parameter container object
    pc = pci;

    // previous functions incremented the indices appropriately, so now we need to get the norm again
    // NOTE: this is the weight for the function about to be called.
    fptype normFactor = pc.getNormalization(0);

    fptype last = callFunction(evt, pc);
    ret += (1 - totalWeight) * last * normFactor;

    return ret;
}

__device__ fptype device_AddPdfsExt(fptype *evt, ParameterContainer &pc) {
    int numParameters  = pc.getNumParameters();
    fptype ret         = 0;
    fptype totalWeight = 0;

    // make a copy of our parameter container
    ParameterContainer pci = pc;

    // We only call increment once we read our weight/norm for the first iteration.
    pci.incrementIndex();

    for(int i = 0; i < numParameters; i++) {
        // grab the weight parameter from addPdf
        fptype weight     = pc.getParameter(i);
        //  Grab the normalization for the specific component
        fptype normFactor = pci.getNormalization(0);

        fptype curr = callFunction(evt, pci);
        ret += weight * curr * normFactor;

        totalWeight += weight;
    }

    pc = pci;
    ret /= totalWeight;

    return ret;
}

__device__ device_function_ptr ptr_to_AddPdfs    = device_AddPdfs;
__device__ device_function_ptr ptr_to_AddPdfsExt = device_AddPdfsExt;

AddPdf::AddPdf(std::string n, std::vector<Variable> weights, std::vector<PdfBase *> comps)
    : CombinePdf("AddPdf", n)
    , extended(true) {
    if(weights.size() != comps.size() && (weights.size() + 1) != comps.size())
        throw GooFit::GeneralError("Size of weights {} (+1) != comps {}", weights.size(), comps.size());

    // Indices stores (function index)(function parameter index)(weight index) triplet for each component.
    // Last component has no weight index unless function is extended.
    for(PdfBase *p : comps) {
        components.push_back(p);
        if(components.back() == nullptr)
            throw GooFit::GeneralError("Invalid component");
    }

    observablesList = getObservables();

    for(unsigned int w = 0; w < weights.size(); ++w) {
        if(components[w] == nullptr)
            throw GooFit::GeneralError("Invalid component");
        registerParameter(weights[w]);
    }

    if(components.back() == nullptr)
        throw GooFit::GeneralError("Invalid component");

    if(weights.size() < components.size()) {
        extended = false;
    }

    if(extended)
        registerFunction("ptr_to_AddPdfsExt", ptr_to_AddPdfsExt);
    else
        registerFunction("ptr_to_AddPdfs", ptr_to_AddPdfs);

    initialize();
}

AddPdf::AddPdf(std::string n, Variable frac1, PdfBase *func1, PdfBase *func2)
    : CombinePdf("AddPdf", n, frac1)
    , extended(false) {
    // Special-case constructor for common case of adding two functions.
    components.push_back(func1);
    components.push_back(func2);

    observablesList = getObservables();

    registerFunction("ptr_to_AddPdfs", ptr_to_AddPdfs);

    initialize();
}

__host__ fptype AddPdf::normalize() {
    // if (cpuDebug & 1) std::cout << "Normalizing AddPdf " << getName() << std::endl;

    fptype ret         = 0;
    fptype totalWeight = 0;

    for(unsigned int i = 0; i < components.size() - 1; ++i) {
        // fptype weight = host_parameters[parametersIdx + 3*i + 1];
        fptype weight = parametersList[i].getValue();
        totalWeight += weight;
        fptype curr = components[i]->normalize();
        ret += curr * weight;
    }

    fptype last = components.back()->normalize();

    if(extended) {
        // fptype lastWeight = host_parameters[parametersIdx + 2];
        fptype lastWeight = parametersList[components.size() - 1];
        totalWeight += lastWeight;
        ret += last * lastWeight;
        ret /= totalWeight;
    } else {
        ret += (1 - totalWeight) * last;
    }

    host_normalizations[normalIdx + 1] = 1.0;
    cachedNormalization                = 1.0;

    if(getCommonNorm()) {
        // Want to normalize this as
        // (f1 A + (1-f1) B) / int (f1 A + (1-f1) B)
        // instead of default
        // (f1 A / int A) + ((1-f1) B / int B).

        for(auto component : components) {
            // host_normalizations[component->getParameterIndex()] = (1.0 / ret);
            component->setNormalization(1.0 / ret);
        }
    }

    // if (cpuDebug & 1) std::cout << getName() << " integral returning " << ret << std::endl;
    return ret;
}

__host__ double AddPdf::calculateNLL() {
    double ret = GooPdf::calculateNLL() / 2.0;

    if(extended) {
        fptype expEvents = 0;

        for(unsigned int i = 0; i < components.size(); ++i) {
            // expEvents += host_parameters[parametersIdx + 3 * (i + 1)];
            expEvents += parametersList[i].getValue();
        }

        // Log-likelihood of numEvents with expectation of exp is (-exp + numEvents*ln(exp) - ln(numEvents!)).
        // The last is constant, so we drop it; and then multiply by minus one to get the negative log-likelihood.
        ret += (expEvents - numEvents * log(expEvents));
    }

    return ret * 2.0;
}
} // namespace GooFit
