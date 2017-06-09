#include "goofit/PDFs/combine/AddPdf.h"
#include "goofit/Error.h"
#include "goofit/detail/ThrustOverride.h"

#include <thrust/iterator/constant_iterator.h>
#include <thrust/transform_reduce.h>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

__device__ fptype device_AddPdfs(fptype *evt, fptype *p, unsigned int *indices) {
    int numParameters  = RO_CACHE(indices[0]);
    fptype ret         = 0;
    fptype totalWeight = 0;

    for(int i = 1; i < numParameters - 3; i += 3) {
        totalWeight += RO_CACHE(p[RO_CACHE(indices[i + 2])]);
        fptype curr   = callFunction(evt, RO_CACHE(indices[i]), RO_CACHE(indices[i + 1]));
        fptype weight = RO_CACHE(p[RO_CACHE(indices[i + 2])]);
        ret += weight * curr * normalisationFactors[RO_CACHE(indices[i + 1])];

        // if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("Add comp %i: %f * %f * %f = %f (%f)\n", i, weight, curr, normalisationFactors[indices[i+1]],
        // weight*curr*normalisationFactors[indices[i+1]], ret);
    }

    // numParameters does not count itself. So the array structure for two functions is
    // nP | F P w | F P
    // in which nP = 5. Therefore the parameter index for the last function pointer is nP, and the function index is
    // nP-1.
    // fptype last = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[numParameters-1]])))(evt, p,
    // paramIndices + indices[numParameters]);
    fptype last = callFunction(evt, RO_CACHE(indices[numParameters - 1]), RO_CACHE(indices[numParameters]));
    ret += (1 - totalWeight) * last * normalisationFactors[RO_CACHE(indices[numParameters])];

    // if ((THREADIDX < 50) && (isnan(ret))) printf("NaN final component %f %f\n", last, totalWeight);

    // if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
    // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    // printf("Add final: %f * %f * %f = %f (%f)\n", (1 - totalWeight), last,
    // normalisationFactors[indices[numParameters]], (1 - totalWeight) *last*
    // normalisationFactors[indices[numParameters]], ret);

    return ret;
}

__device__ fptype device_AddPdfsExt(fptype *evt, fptype *p, unsigned int *indices) {
    // numParameters does not count itself. So the array structure for two functions is
    // nP | F P w | F P w
    // in which nP = 6.

    int numParameters  = RO_CACHE(indices[0]);
    fptype ret         = 0;
    fptype totalWeight = 0;

    for(int i = 1; i < numParameters; i += 3) {
        // fptype curr = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[i]])))(evt, p,
        // paramIndices + indices[i+1]);
        fptype curr   = callFunction(evt, RO_CACHE(indices[i]), RO_CACHE(indices[i + 1]));
        fptype weight = RO_CACHE(p[RO_CACHE(indices[i + 2])]);
        ret += weight * curr * normalisationFactors[RO_CACHE(indices[i + 1])];

        totalWeight += weight;
        // if ((gpuDebug & 1) && (THREADIDX == 0) && (0 == BLOCKIDX))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("AddExt: %i %E %f %f %f %f %f %f\n", i, curr, weight, ret, totalWeight,
        // normalisationFactors[indices[i+1]], evt[0], evt[8]);
    }

    ret /= totalWeight;
    // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    // if ((gpuDebug & 1) && (THREADIDX == 0) && (0 == BLOCKIDX))
    // printf("AddExt result: %f\n", ret);

    return ret;
}

__device__ device_function_ptr ptr_to_AddPdfs    = device_AddPdfs;
__device__ device_function_ptr ptr_to_AddPdfsExt = device_AddPdfsExt;

AddPdf::AddPdf(std::string n, std::vector<Variable *> weights, std::vector<PdfBase *> comps)
    : GooPdf(nullptr, n)
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

    observables = getObservables();

    std::vector<unsigned int> pindices;

    for(unsigned int w = 0; w < weights.size(); ++w) {
        if(components[w] == nullptr)
            throw GooFit::GeneralError("Invalid component");
        pindices.push_back(components[w]->getFunctionIndex());
        pindices.push_back(components[w]->getParameterIndex());
        pindices.push_back(registerParameter(weights[w]));
    }

    if(components.back() == nullptr)
        throw GooFit::GeneralError("Invalid component");

    if(weights.size() < components.size()) {
        pindices.push_back(components.back()->getFunctionIndex());
        pindices.push_back(components.back()->getParameterIndex());
        extended = false;
    }

    if(extended)
        GET_FUNCTION_ADDR(ptr_to_AddPdfsExt);
    else
        GET_FUNCTION_ADDR(ptr_to_AddPdfs);

    initialize(pindices);
}

AddPdf::AddPdf(std::string n, Variable *frac1, PdfBase *func1, PdfBase *func2)
    : GooPdf(nullptr, n)
    , extended(false) {
    // Special-case constructor for common case of adding two functions.
    components.push_back(func1);
    components.push_back(func2);
    observables = getObservables();

    std::vector<unsigned int> pindices;
    pindices.push_back(func1->getFunctionIndex());
    pindices.push_back(func1->getParameterIndex());
    pindices.push_back(registerParameter(frac1));

    pindices.push_back(func2->getFunctionIndex());
    pindices.push_back(func2->getParameterIndex());

    GET_FUNCTION_ADDR(ptr_to_AddPdfs);

    initialize(pindices);
}

__host__ fptype AddPdf::normalize() const {
    // if (cpuDebug & 1) std::cout << "Normalising AddPdf " << getName() << std::endl;

    fptype ret         = 0;
    fptype totalWeight = 0;

    for(unsigned int i = 0; i < components.size() - 1; ++i) {
        fptype weight = host_params[host_indices[parameters + 3 * (i + 1)]];
        totalWeight += weight;
        fptype curr = components[i]->normalize();
        ret += curr * weight;
    }

    fptype last = components.back()->normalize();

    if(extended) {
        fptype lastWeight = host_params[host_indices[parameters + 3 * components.size()]];
        totalWeight += lastWeight;
        ret += last * lastWeight;
        ret /= totalWeight;
    } else {
        ret += (1 - totalWeight) * last;
    }

    host_normalisation[parameters] = 1.0;

    if(getSpecialMask() & PdfBase::ForceCommonNorm) {
        // Want to normalize this as
        // (f1 A + (1-f1) B) / int (f1 A + (1-f1) B)
        // instead of default
        // (f1 A / int A) + ((1-f1) B / int B).

        for(auto component : components) {
            host_normalisation[component->getParameterIndex()] = (1.0 / ret);
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
            expEvents += host_params[host_indices[parameters + 3 * (i + 1)]];
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
