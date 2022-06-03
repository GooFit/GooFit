#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/utilities/VariableBinTransform1DPdf.h>

namespace GooFit {

__device__ auto device_VarBinTransform1D(fptype *evt, ParameterContainer &pc) -> fptype {
    // Index structure: nP lim1 lim2 ...
    int ret = 0;
    // int previousSize = 1;
    // printf("[%i, %i] Bin Transform: %i %i %f %f\n", THREADIDX, BLOCKIDX, numObservables, previousSize, evt[0],
    // evt[1]);
    int id          = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    fptype obsValue = evt[id];
    if(obsValue < 0)
        obsValue = -obsValue;
    int numLimits = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    for(int i = 0; i < numLimits; ++i) {
        fptype lowerLimit = RO_CACHE(pc.constants[pc.constantIdx + i + 2]);
        if(obsValue < lowerLimit)
            break;
        ret++;
    }

    return fptype(ret);
}

__device__ device_function_ptr ptr_to_VarBinTransform1D = device_VarBinTransform1D;

// Notice that bin sizes and limits can be different, for this purpose, than what's implied by the Variable members.
__host__
VariableBinTransform1DPdf::VariableBinTransform1DPdf(std::string n, Observable _x, std::vector<fptype> binlimits)
    : GooPdf("VariableBinTransform1DPdf", n, _x) {
    unsigned int numLimits = binlimits.size(); // Excluding the min & max values for _x

    constantsList.push_back(numLimits);
    for(size_t i = 0; i < numLimits; ++i) {
        constantsList.push_back(binlimits[i]);
    }

    registerFunction("ptr_to_VarBinTransform1D", ptr_to_VarBinTransform1D);
    // initialize() ?
}

} // namespace GooFit
