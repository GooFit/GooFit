#include <goofit/PDFs/utility/VariableBinTransform1DPdf.h>

using namespace std;

namespace GooFit {

__device__ fptype device_VarBinTransform1D(fptype *evt, fptype *p, unsigned int *indices) {
    // Index structure: nP lim1 lim2 ...
    int ret = 0;
    // int previousSize = 1;
    // printf("[%i, %i] Bin Transform: %i %i %f %f\n", THREADIDX, BLOCKIDX, numObservables, previousSize, evt[0],
    // evt[1]);
    fptype obsValue = evt[indices[2 + indices[0]]];
    if(obsValue < 0)
        obsValue = -obsValue;
    int numLimits = indices[1];
    for(int i = 0; i < numLimits; ++i) {
        fptype lowerLimit = functorConstants[indices[i + 2]];
        if(obsValue < lowerLimit)
            break;
        ret++;
    }

    return fptype(ret);
}

__device__ device_function_ptr ptr_to_VarBinTransform1D = device_VarBinTransform1D;

// Notice that bin sizes and limits can be different, for this purpose, than what's implied by the Variable members.
__host__ VariableBinTransform1DPdf::VariableBinTransform1DPdf(std::string n, Observable _x, vector<fptype> binlimits)
    : GooPdf(n, _x) {
    unsigned int numLimits = binlimits.size(); // Excluding the min & max values for _x
    cIndex                 = registerConstants(numLimits);
    std::vector<fptype> host_constants;
    std::vector<unsigned int> pindices{numLimits};
    for(size_t i = 0; i < numLimits; ++i) {
        pindices.push_back(cIndex + i);
        host_constants.push_back(binlimits[i]); // cIndex will be accounted for by offset in memcpy
    }

    MEMCPY_TO_SYMBOL(functorConstants,
                     host_constants.data(),
                     numLimits * sizeof(fptype),
                     cIndex * sizeof(fptype),
                     cudaMemcpyHostToDevice);

    GET_FUNCTION_ADDR(ptr_to_VarBinTransform1D);
    initialize(pindices);
}

} // namespace GooFit
