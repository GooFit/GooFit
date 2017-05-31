#include "goofit/PDFs/combine/MappedPdf.h"

namespace GooFit {

__device__ fptype device_Mapped(fptype *evt, fptype *p, unsigned int *indices) {
    // Structure : nP mapFunctionIndex mapParamIndex functionIndex1 parameterIndex1 functionIndex2 parameterIndex2 ...

    // Find mapping between event variables and function to evaluate
    unsigned int mapFunction = RO_CACHE(indices[1]);
    // This is an index into the MappedPdf's list of functions
    // int targetFunction = (int) floor(0.5 +
    // (*(reinterpret_cast<device_function_ptr>(device_function_table[mapFunction])))(evt, p, paramIndices +
    // indices[2]));
    auto targetFunction = static_cast<int>(floor(0.5 + callFunction(evt, mapFunction, RO_CACHE(indices[2]))));

    targetFunction *= 2; // Because there are two pieces of information about each function
    targetFunction += 3; // Because first function information begins at index 3

    // fptype ret = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[targetFunction]])))(evt, p,
    // paramIndices + indices[targetFunction + 1]);
    fptype ret = callFunction(evt, RO_CACHE(indices[targetFunction]), RO_CACHE(indices[targetFunction + 1]));
    ret *= normalisationFactors[RO_CACHE(indices[targetFunction + 1])];
    // if (gpuDebug & 1)
    // if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX))
    // printf("[%i, %i] Mapped: %i (%f %f %f %f) %f\n", BLOCKIDX, THREADIDX, targetFunction, evt[0], evt[1], evt[2],
    // evt[3], ret);
    return ret;
}

__device__ device_function_ptr ptr_to_Mapped = device_Mapped;

__host__ MappedPdf::MappedPdf(std::string n, GooPdf *m, std::vector<GooPdf *> &t)
    : GooPdf(nullptr, n) {
    components.push_back(m);
    std::vector<unsigned int> pindices;
    pindices.push_back(m->getFunctionIndex());
    pindices.push_back(m->getParameterIndex());

    std::set<int> functionIndicesUsed;

    for(GooPdf *f : t) {
        components.push_back(f);
        pindices.push_back(f->getFunctionIndex());
        pindices.push_back(f->getParameterIndex());
        functionIndicesUsed.insert(f->getFunctionIndex());
    }

    if(functionIndicesUsed.size() > 1) {
        std::cout << "Warning: More than one function type given to MappedPdf " << getName()
                  << " constructor. This may slow execution by causing sequential evaluations.\n";
    }

    observables = getObservables();
    GET_FUNCTION_ADDR(ptr_to_Mapped);
    initialize(pindices);
}

__host__ fptype MappedPdf::normalize() const {
    // std::cout << "Normalising MappedPdf " << getName() << std::endl;
    fptype ret = 0;

    for(unsigned int i = 1; i < components.size(); ++i) { // No need to normalize mapping function.
        fptype curr = components[i]->normalize();
        ret += curr;
    }

    host_normalisation[parameters] = 1.0;
    return ret;
}
} // namespace GooFit
