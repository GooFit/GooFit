#include "goofit/PDFs/basic/BinTransformPdf.h"

namespace GooFit {

__device__ fptype device_BinTransform(fptype *evt, ParameterContainer &pc) {
    // Index structure: nP lim1 bin1 lim2 bin2 ... nO o1 o2
    int numConstants = pc.constants[pc.constantIdx];
    int numObservables = pc.observables[pc.observableIdx];
    int ret            = 0;
    int previousSize   = 1;

    for(int i = 0; i < numObservables; ++i) {
        int id = pc.observables[pc.observableIdx + i + 1];
        fptype obsValue   = evt[id];
        fptype lowerLimit = pc.constants[pc.constantIdx + i*3 + 1];
        fptype binSize    = pc.constants[pc.constantIdx + i*3 + 2];
        int numBins       = pc.constants[pc.constantIdx + i*3 + 3];

        auto localBin = static_cast<int>(floor((obsValue - lowerLimit) / binSize));
        ret += localBin * previousSize;
        previousSize *= numBins;
    }

    pc.incrementIndex(1, 0, numConstants, numObservables, 1);

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
    //cIndex               = registerConstants(2 * obses.size());
    //auto *host_constants = new fptype[2 * obses.size()];
    std::vector<unsigned int> pindices;

    //setup the observables
    for(unsigned int i = 0; i < obses.size(); ++i)
        registerObservable(obses[i]);

    observablesList = getObservables ();

    //add limits for each observable (variable)
    for(unsigned int i = 0; i < obses.size(); ++i) {
        constantsList.push_back (limits[i]);
        constantsList.push_back (binSizes[i]);
        constantsList.push_back (numBins[i]);
    }

    GET_FUNCTION_ADDR(ptr_to_BinTransform);
    initialize(pindices);
}

__host__ void BinTransformPdf::recursiveSetIndices () {
    GET_FUNCTION_ADDR(ptr_to_BinTransform);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_BinTransform");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;

    populateArrays ();
}

} // namespace GooFit
