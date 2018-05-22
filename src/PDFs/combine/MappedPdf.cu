#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/combine/MappedPdf.h>

namespace GooFit {

__device__ fptype device_Mapped(fptype *evt, ParameterContainer &pc) {
    // Structure : nP mapFunctionIndex mapParamIndex functionIndex1 parameterIndex1 functionIndex2 parameterIndex2 ...

    // Find mapping between event variables and function to evaluate
    unsigned int numTargets = pc.getConstant(0);

    // Mapping PDF happens directly after, so just increment.
    // pc.incrementIndex (1, 0, 1, 0, 1);
    pc.incrementIndex();

    // This is an index into the MappedPdf's list of functions
    // int targetFunction = (int) floor(0.5 +
    // (*(reinterpret_cast<device_function_ptr>(device_function_table[mapFunction])))(evt, p, paramIndices +
    // indices[2]));
    auto targetFunction = static_cast<int>(floor(0.5 + callFunction(evt, pc)));

    // targetFunction *= 2; // Because there are two pieces of information about each function
    // targetFunction += 3; // Because first function information begins at index 3

    // Unsure these carry over...
    // targetFunction *= 2;
    // targetFunction += 3;

    // numTargets; increment past our set of mapping functions and our target functions after it is handled.
    unsigned int funcIdx = pc.funcIdx;

    while(pc.funcIdx < funcIdx + targetFunction)
        pc.incrementIndex();

    // fptype ret = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[targetFunction]])))(evt, p,
    // paramIndices + indices[targetFunction + 1]);
    fptype norm = pc.getNormalisation(0);
    fptype ret  = callFunction(evt, pc);
    ret *= norm;

    // increment our functions here...
    while(pc.funcIdx < funcIdx + numTargets)
        pc.incrementIndex();

    // if (gpuDebug & 1)
    // if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX))
    // printf("[%i, %i] Mapped: %i (%f %f %f %f) %f\n", BLOCKIDX, THREADIDX, targetFunction, evt[0], evt[1], evt[2],
    // evt[3], ret);
    return ret;
}

__device__ device_function_ptr ptr_to_Mapped = device_Mapped;

__host__ MappedPdf::MappedPdf(std::string n, GooPdf *m, std::vector<GooPdf *> &t)
    : GooPdf(n) {
    components.push_back(m);

    std::set<int> functionIndicesUsed;

    for(GooPdf *f : t) {
        components.push_back(f);
        // pindices.push_back(f->getFunctionIndex());
        // pindices.push_back(f->getParameterIndex());
        // functionIndicesUsed.insert(f->getFunctionIndex());
    }

    // if(functionIndicesUsed.size() > 1) {
    //    std::cout << "Warning: More than one function type given to MappedPdf " << getName()
    //              << " constructor. This may slow execution by causing sequential evaluations.\n";
    //}

    // This makes sure we have the appropriate amount of obs in our structure
    observablesList = getObservables();

    // add a constant value for the number of 't' functions, skipping 'm'.
    registerConstant(components.size() - 1);

    registerFunction("ptr_to_Mapped", ptr_to_Mapped);

    initialize();
}

__host__ fptype MappedPdf::normalize() const {
    // std::cout << "Normalising MappedPdf " << getName() << std::endl;
    fptype ret = 0;

    for(unsigned int i = 1; i < components.size(); ++i) { // No need to normalize mapping function.
        fptype curr = components[i]->normalize();
        ret += curr;
    }

    host_normalisations[normalIdx + 1] = 1.0;
    return ret;
}
} // namespace GooFit
