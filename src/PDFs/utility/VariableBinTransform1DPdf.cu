#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/utility/VariableBinTransform1DPdf.h>

using namespace std;

namespace GooFit {

__device__ fptype device_VarBinTransform1D(fptype *evt, ParameterContainer &pc) {
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
__host__ VariableBinTransform1DPdf::VariableBinTransform1DPdf(std::string n, Observable _x, vector<fptype> binlimits)
    : GooPdf(n, _x) {
    unsigned int numLimits = binlimits.size(); // Excluding the min & max values for _x
    // cIndex                 = registerConstants(numLimits);
    // std::vector<fptype> host_constants;
    // std::vector<unsigned int> pindices{numLimits};
    // for(size_t i = 0; i < numLimits; ++i) {
    //    pindices.push_back(cIndex + i);
    //    host_constants.push_back(binlimits[i]); // cIndex will be accounted for by offset in memcpy
    //}

    constantsList.push_back(numLimits);
    for(size_t i = 0; i < numLimits; ++i) {
        constantsList.push_back(binlimits[i]);
    }

    // MEMCPY_TO_SYMBOL(functorConstants,
    //                 host_constants.data(),
    //                 numLimits * sizeof(fptype),
    //                 cIndex * sizeof(fptype),
    //                 cudaMemcpyHostToDevice);

    // GET_FUNCTION_ADDR(ptr_to_VarBinTransform1D);
    // initialize(pindices);
}

void VariableBinTransform1DPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_VarBinTransform1D);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_VarBinTransform1D");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

} // namespace GooFit
