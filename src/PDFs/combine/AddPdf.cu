#include <goofit/Error.h>
#include <goofit/Log.h>
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
        fptype norm = pci.getNormalisation(0);

        // call the first function to add in our PDF.
        fptype curr = callFunction(evt, pci);

        ret += weight * curr * norm;
    }

    // restore our new parameter container object
    pc = pci;

    // previous functions incremented the indices appropriately, so now we need to get the norm again
    // NOTE: this is the weight for the function about to be called.
    fptype normFactor = pc.getNormalisation(0);

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
        // grab the weight value
        fptype weight     = pci.getParameter(i);
        fptype normFactor = pci.getNormalisation(0);

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
    : GooPdf(n)
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

    initialize();
}

AddPdf::AddPdf(std::string n, Variable frac1, PdfBase *func1, PdfBase *func2)
    : GooPdf(n)
    , extended(false) {
    // Special-case constructor for common case of adding two functions.
    components.push_back(func1);
    components.push_back(func2);

    observablesList = getObservables();

    registerParameter(frac1);

    initialize();
}

__host__ void AddPdf::recursiveSetIndices() {
    if(extended) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_AddPdfsExt");
        GET_FUNCTION_ADDR(ptr_to_AddPdfsExt);
    } else {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_AddPdfs");
        GET_FUNCTION_ADDR(ptr_to_AddPdfs);
    }

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

__host__ fptype AddPdf::normalize() const {
    // if (cpuDebug & 1) std::cout << "Normalising AddPdf " << getName() << std::endl;

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
        fptype lastWeight = host_parameters[parametersIdx + 2];
        totalWeight += lastWeight;
        ret += last * lastWeight;
        ret /= totalWeight;
    } else {
        ret += (1 - totalWeight) * last;
    }

    host_normalisations[normalIdx + 1] = 1.0;

    // TODO: Unsure of the exact location for this normalise...
    if(getSpecialMask() & PdfBase::ForceCommonNorm) {
        // Want to normalize this as
        // (f1 A + (1-f1) B) / int (f1 A + (1-f1) B)
        // instead of default
        // (f1 A / int A) + ((1-f1) B / int B).

        for(unsigned int i = 0; i < components.size(); ++i) {
            host_normalisations[components[i]->getParameterIndex()] = (1.0 / ret);
        }
    }

    // if (cpuDebug & 1) std::cout << getName() << " integral returning " << ret << std::endl;
    return ret;
}

__host__ double AddPdf::sumOfNll(int numVars) const {
    static thrust::plus<double> cudaPlus;
    thrust::constant_iterator<int> eventSize(numVars);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    double dummy = 0;

    thrust::counting_iterator<int> eventIndex(0);

    double ret;
#ifdef GOOFIT_MPI
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    goofit_policy my_policy;
    double r = thrust::transform_reduce(
        my_policy,
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#else
    double r = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#endif

    MPI_Allreduce(&r, &ret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    goofit_policy my_policy;
    ret = thrust::transform_reduce(
        my_policy,
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#else
    ret = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#endif
#endif

    if(extended) {
        fptype expEvents = 0;

        // std::cout << "Weights:";
        for(unsigned int i = 0; i < components.size(); ++i) {
            expEvents += host_parameters[parametersIdx + 3 * (i + 1)];
            // std::cout << " " << host_params[host_indices[parameters + 3*(i+1)]];
        }

        // Log-likelihood of numEvents with expectation of exp is (-exp + numEvents*ln(exp) - ln(numEvents!)).
        // The last is constant, so we drop it; and then multiply by minus one to get the negative log-likelihood.
        ret += (expEvents - numEvents * log(expEvents));
        // std::cout << " " << expEvents << " " << numEvents << " " << (expEvents - numEvents*log(expEvents)) <<
        // std::endl;
    }

    return ret;
}
} // namespace GooFit
