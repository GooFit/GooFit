#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/BinTransformPdf.h>

namespace GooFit {

__device__ auto device_BinTransform(fptype *evt, ParameterContainer &pc) -> fptype {
    // Index structure: nP lim1 bin1 lim2 bin2 ... nO o1 o2
    int numConstants   = pc.getNumConstants();
    int numObservables = pc.getNumObservables();
    int ret            = 0;
    int previousSize   = 1;

    for(int i = 0; i < numObservables; ++i) {
        int id            = pc.getObservable(i);
        fptype obsValue   = RO_CACHE(evt[id]);
        fptype lowerLimit = pc.getConstant(i * 3);
        fptype binSize    = pc.getConstant(i * 3 + 1);
        int numBins       = pc.getConstant(i * 3 + 2);

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
                                          std::vector<Observable> obses,
                                          std::vector<fptype> limits,
                                          std::vector<fptype> binSizes,
                                          std::vector<int> numBins)
    : GooPdf("BinTransformPdf", n) {
    // setup the observables
    for(const auto &obse : obses)
        registerObservable(obse);

    observablesList = getObservables();

    // add limits for each observable (variable)
    for(unsigned int i = 0; i < obses.size(); ++i) {
        registerConstant(limits[i]);
        registerConstant(binSizes[i]);
        registerConstant(numBins[i]);
    }

    registerFunction("ptr_to_BinTransform", ptr_to_BinTransform);

    initialize();
}

} // namespace GooFit
