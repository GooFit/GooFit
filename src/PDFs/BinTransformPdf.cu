#include "goofit/PDFs/basic/BinTransformPdf.h"

namespace GooFit {

__device__ fptype device_BinTransform(fptype *evt, fptype *p, unsigned int *indices) {
    // Index structure: nP lim1 bin1 lim2 bin2 ... nO o1 o2
    int numObservables = indices[1 + indices[0]];
    int ret            = 0;
    int previousSize   = 1;

    // printf("[%i, %i] Bin Transform: %i %i %f %f\n", THREADIDX, BLOCKIDX, numObservables, previousSize, evt[0],
    // evt[1]);
    for(int i = 0; i < numObservables; ++i) {
        fptype obsValue   = evt[indices[2 + indices[0] + i]];
        fptype lowerLimit = functorConstants[indices[i * 3 + 1]];
        fptype binSize    = functorConstants[indices[i * 3 + 2]];
        int numBins       = indices[i * 3 + 3];

        auto localBin = static_cast<int>(floor((obsValue - lowerLimit) / binSize));
        ret += localBin * previousSize;
        previousSize *= numBins;
    }

    return fptype(ret);
}

__device__ device_function_ptr ptr_to_BinTransform = device_BinTransform;

// Notice that bin sizes and limits can be different, for this purpose, than what's implied by the Variable members.
__host__ BinTransformPdf::BinTransformPdf(std::string n,
                                          std::vector<Variable *> obses,
                                          std::vector<fptype> limits,
                                          std::vector<fptype> binSizes,
                                          std::vector<int> numBins)
    : GooPdf(nullptr, n) {
    cIndex               = registerConstants(2 * obses.size());
    auto *host_constants = new fptype[2 * obses.size()];
    std::vector<unsigned int> pindices;

    for(unsigned int i = 0; i < obses.size(); ++i) {
        registerObservable(obses[i]);
        pindices.push_back(cIndex + 2 * i);
        pindices.push_back(cIndex + 2 * i + 1);
        pindices.push_back(numBins[i]);

        host_constants[2 * i]     = limits[i]; // cIndex will be accounted for by offset in memcpy
        host_constants[2 * i + 1] = binSizes[i];
    }

    MEMCPY_TO_SYMBOL(functorConstants,
                     host_constants,
                     2 * obses.size() * sizeof(fptype),
                     cIndex * sizeof(fptype),
                     cudaMemcpyHostToDevice);
    delete[] host_constants;

    GET_FUNCTION_ADDR(ptr_to_BinTransform);
    initialize(pindices);
}

} // namespace GooFit
