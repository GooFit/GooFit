#include <algorithm>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/InterHistPdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__constant__ fptype *dev_base_interhists[100]; // Multiple histograms for the case of multiple PDFs
#define OBS_CODE 4242424242
// This number is presumably so high that it will never collide
// with an actual parameter index. It indicates that this dimension
// is an event observable.

// dev_powi is implemented in SmoothHistogramPdf.cu.

__device__ auto device_InterHistogram(fptype *evt, ParameterContainer &pc) -> fptype {
    // Structure is
    // nP totalHistograms (idx1 limit1 step1 bins1) (idx2 limit2 step2 bins2) nO o1 o2
    // where limit and step are indices into functorConstants.

    int numVars          = int(pc.getConstant(1) - 1) / 4; //(indices[0] - 1) / 4;
    int globalBin        = 0;
    int previous         = 1;
    int myHistogramIndex = pc.getConstant(0); // indices[1];
    fptype binDistances[10];                  // Ten dimensions should be more than enough!
    // Distance from bin center in units of bin width in each dimension.

    for(int i = 0; i < numVars; ++i) {
        fptype currVariable   = 0;
        unsigned int varIndex = pc.getConstant(2 + i * 4); // constantindices[2 + 4 * i];

        // check where we get our value
        if(varIndex == OBS_CODE) {
            // Interpret this number as observable index.
            // Notice that this if does not cause a fork
            // - all threads will hit the same index and
            // make the same decision.
            int id       = pc.getObservable(i);
            currVariable = RO_CACHE(evt[id]); // evt[indices[indices[0] + 2 + observablesSeen++]];
        } else {
            // Interpret as parameter index.
            currVariable = pc.getParameter(varIndex);
        }

        int lowerBoundIdx = 3 + 4 * i;
        fptype lowerBound = pc.getParameter(3 + i * 4);     // functorConstants[indices[lowerBoundIdx + 0]];
        fptype step       = pc.getParameter(4 + i * 4 + 2); // functorConstants[indices[lowerBoundIdx + 1]];

        currVariable -= lowerBound;
        currVariable /= step;

        auto localBin   = static_cast<int>(floor(currVariable));
        binDistances[i] = currVariable - localBin - fptype(0.5);
        globalBin += previous * localBin;
        previous *= pc.getConstant(lowerBoundIdx + 1); // indices[lowerBoundIdx + 2];

        if(0 == THREADIDX + BLOCKIDX)
            printf("Variable %i: %f %f %i\n", i, currVariable, currVariable * step + lowerBound, localBin);
    }

    fptype *myHistogram = dev_base_interhists[myHistogramIndex];
    fptype ret          = 0;

    //------------------
    //     |     |     |
    //  3  |  4  |  5  |
    //     |     |     |
    //------------------
    //    x|     |     |
    //  0  |  1  |  2  |
    //     |     |     |
    //------------------

    fptype totalWeight = 0;
    int totalBins      = dev_powi(3, numVars);

    for(int i = 0; i < totalBins; ++i) {
        int currBin          = globalBin;
        int localPrevious    = 1;
        int trackingBin      = globalBin;
        bool offSomeAxis     = false;
        fptype currentWeight = 0;

        // Loop over vars to get offset for each one.
        for(int v = 0; v < numVars; ++v) {
            int localNumBins = pc.parameters[pc.parameterIdx + 1 + v]; // indices[4 * (v + 1) + 1];
            int offset       = ((i / dev_powi(3, v)) % 3) - 1;

            currBin += offset * localPrevious;
            localPrevious *= localNumBins;

            int currVarBin = trackingBin % localNumBins;
            trackingBin /= localNumBins;

            if(currVarBin + offset < 0)
                offSomeAxis = true;

            if(currVarBin + offset >= localNumBins)
                offSomeAxis = true;

            fptype currDist = binDistances[v];
            currDist -= offset;
            currentWeight += currDist * currDist;

            if(0 == THREADIDX + BLOCKIDX)
                printf("%i, %i: %f %f %f %i %s\n",
                       i,
                       v,
                       currDist,
                       binDistances[v],
                       currentWeight,
                       offset,
                       offSomeAxis ? "off" : "on");
        }

        // Only interpolate the four closest boxes (in two dimensions; more in three dimensions).
        currentWeight       = currentWeight > 0
                                  ? (currentWeight <= sqrt(static_cast<fptype>(numVars)) ? 1 / sqrt(currentWeight) : 0)
                                  : 0;
        fptype currentEntry = offSomeAxis ? 0 : myHistogram[currBin];
        ret += currentWeight * currentEntry;
        totalWeight += currentWeight;

        if(0 == THREADIDX + BLOCKIDX)
            printf(
                "Adding bin content %i %f with weight %f for total %f.\n", currBin, currentEntry, currentWeight, ret);
    }

    ret /= totalWeight;
    return ret;
}

__device__ device_function_ptr ptr_to_InterHistogram = device_InterHistogram;

__host__
InterHistPdf::InterHistPdf(std::string n, BinnedDataSet *x, std::vector<Variable> params, std::vector<Observable> obses)
    : GooPdf("InterHistPdf", n)
    , numVars(x->numVariables()) {
    static unsigned int totalHistograms = 0;
    totalEvents                         = 0;

    // push on the histogram index and number of variables.
    registerConstant(totalHistograms);
    registerConstant(numVars);

    int varIndex = 0;

    for(Observable var : x->getObservables()) {
        registerObservable(var);
        registerConstant(OBS_CODE);

        registerConstant(var.getLowerLimit());
        registerConstant(var.getBinSize());
        registerConstant(var.getNumBins());

        // NB, do not put cIndex here, it is accounted for by the offset in MEMCPY_TO_SYMBOL below.
        varIndex++;
    }

    unsigned int numbins = x->getNumBins();
    thrust::host_vector<fptype> host_histogram;

    for(unsigned int i = 0; i < numbins; ++i) {
        fptype curr = x->getBinContent(i);
        host_histogram.push_back(curr);
        totalEvents += curr;
    }

    dev_base_histogram = new thrust::device_vector<fptype>(host_histogram);
    static fptype *dev_address[1];
    dev_address[0] = (&((*dev_base_histogram)[0])).get();
    MEMCPY_TO_SYMBOL(
        dev_base_interhists, dev_address, sizeof(fptype *), totalHistograms * sizeof(fptype *), cudaMemcpyHostToDevice);

    registerFunction("ptr_to_InterHistogram", ptr_to_InterHistogram);

    initialize();

    totalHistograms++;
}

} // namespace GooFit
