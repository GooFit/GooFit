#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__constant__ fptype *dev_base_histograms[100]; // Multiple histograms for the case of multiple PDFs
__constant__ fptype *dev_smoothed_histograms[100];
unsigned int SmoothHistogramPdf::totalHistograms = 0;

__device__ auto device_EvalHistogram(fptype *evt, ParameterContainer &pc) -> fptype {
    // Structure is
    // nP smoothingIndex totalHistograms (limit1 step1 bins1) (limit2 step2 bins2) nO o1 o2
    // where limit and step are indices into functorConstants.

    int numCons          = pc.getNumConstants();
    int numObs           = pc.getNumObservables();
    int numParms         = pc.getNumParameters();
    int numVars          = pc.getConstant(0);
    int globalBinNumber  = 0;
    int previous         = 1;
    int myHistogramIndex = pc.getConstant(1); // 1 only used for smoothing

    for(int i = 0; i < numVars; ++i) {
        int varIndex      = pc.getObservable(i);
        int lowerBoundIdx = 3 * (i + 1);

        fptype currVariable = RO_CACHE(evt[varIndex]);
        fptype lowerBound   = pc.getConstant(i * 3 + 4);
        fptype step         = pc.getConstant(i * 3 + 5);

        currVariable -= lowerBound;
        currVariable /= step;

        auto localBinNumber = static_cast<int>(floor(currVariable));
        globalBinNumber += previous * localBinNumber;

        // This is accessing too far ahead?
        int offset = pc.getConstant(lowerBoundIdx);
        previous *= offset;
    }

    fptype *myHistogram = dev_smoothed_histograms[myHistogramIndex];
    fptype ret          = myHistogram[globalBinNumber];

    pc.incrementIndex(1, numParms, numCons, numObs, 1);

    return ret;
}

struct Smoother {
    int funcIdx;

    __device__ auto operator()(int globalBin) -> fptype {
        ParameterContainer pc;

        while(pc.funcIdx < funcIdx)
            pc.incrementIndex();

        int numVars         = pc.getConstant(0);
        fptype smoothing    = pc.getParameter(0);
        int histIndex       = pc.getConstant(1);
        fptype *myHistogram = dev_base_histograms[histIndex];
        fptype centralValue = myHistogram[globalBin];

        fptype otherBinsTotal  = 0;
        int numSurroundingBins = 0;
        int otherBins          = dev_powi(3, numVars);

        for(int i = 0; i < otherBins; ++i) {
            int currBin       = globalBin;
            int localPrevious = 1;
            int trackingBin   = globalBin;
            bool offSomeAxis  = false;

            for(int v = 0; v < numVars; ++v) {
                int localNumBins = pc.getConstant(3 * (v + 1));
                int offset       = ((i / dev_powi(3, v)) % 3) - 1;

                currBin += offset * localPrevious;
                localPrevious *= localNumBins;

                int currVarBin = trackingBin % localNumBins;
                trackingBin /= localNumBins;

                if(currVarBin + offset < 0)
                    offSomeAxis = true;

                if(currVarBin + offset >= localNumBins)
                    offSomeAxis = true;
            }

            if(currBin == globalBin)
                continue;

            if(offSomeAxis)
                continue; // Out of bounds

            numSurroundingBins++;

            otherBinsTotal += myHistogram[currBin];
        }

        centralValue += otherBinsTotal * smoothing;
        centralValue /= (1 + numSurroundingBins * smoothing);
        return centralValue;
    }
};

__device__ device_function_ptr ptr_to_EvalHistogram = device_EvalHistogram;

__host__ SmoothHistogramPdf::SmoothHistogramPdf(std::string n, BinnedDataSet *hist, Variable smoothing)
    : GooPdf("SmoothHistogramPdf", n, smoothing) {
    int numVars = hist->numVariables();
    totalEvents = 0;

    int varIndex = 0;

    registerConstant(hist->getObservables().size());
    registerConstant(totalHistograms);
    registerConstant(numVars);

    for(const Observable &var : hist->getObservables()) {
        registerObservable(var);
        registerConstant(var.getNumBins());
        registerConstant(var.getLowerLimit());
        registerConstant(var.getBinSize());

        // host_constants[2 * varIndex + 0] = var->getLowerLimit(); // NB, do not put cIndex here, it is accounted for
        // by
        // the offset in MEMCPY_TO_SYMBOL below.
        varIndex++;
    }

    unsigned int numbins = hist->getNumBins();
    thrust::host_vector<fptype> host_histogram;

    for(unsigned int i = 0; i < numbins; ++i) {
        fptype curr = hist->getBinContent(i);
        host_histogram.push_back(curr);
        totalEvents += curr;
    }

    if(totalEvents > 0)
        copyHistogramToDevice(host_histogram);
    else
        std::cout << "Warning: Empty histogram supplied to " << getName()
                  << " not copied to device. Expect copyHistogramToDevice call later.\n";

    registerFunction("ptr_to_EvalHistogram", ptr_to_EvalHistogram);

    initialize();
}

auto pointerToFirst(thrust::device_vector<fptype> *hist) -> fptype * { return (&((*hist)[0])).get(); }

auto pointerToFirst(thrust::host_vector<fptype> *hist) -> fptype * {
    // (*hist) is the host_vector.
    // (*hist)[0] is a 'reference' - Thrust class, not ordinary C++ reference -
    // to the first element of the vector.
    // &((*hist)[0]) is a 'Pointer', as defined by the host_vector, to the location
    // of the 'reference'. Fortunately this is by default fptype*!
    return &((*hist)[0]);
}

__host__ void SmoothHistogramPdf::copyHistogramToDevice(thrust::host_vector<fptype> &host_histogram) {
    dev_base_histogram     = new thrust::device_vector<fptype>(host_histogram);
    dev_smoothed_histogram = new thrust::device_vector<fptype>(host_histogram);
    static fptype *dev_address[1];
    dev_address[0] = pointerToFirst(dev_base_histogram);
    MEMCPY_TO_SYMBOL(
        dev_base_histograms, dev_address, sizeof(fptype *), totalHistograms * sizeof(fptype *), cudaMemcpyHostToDevice);
    dev_address[0] = pointerToFirst(dev_smoothed_histogram);
    MEMCPY_TO_SYMBOL(dev_smoothed_histograms,
                     dev_address,
                     sizeof(fptype *),
                     totalHistograms * sizeof(fptype *),
                     cudaMemcpyHostToDevice);

    totalHistograms++;

    int expectedBins = 1;

    for(Observable &observable : observablesList) {
        expectedBins *= observable.getNumBins();
    }

    if(expectedBins != host_histogram.size()) {
        std::cout << "Warning: Histogram supplied to " << getName() << " has " << host_histogram.size()
                  << " bins, expected " << expectedBins << " - may indicate a problem.\n";
    }
}

__host__ auto SmoothHistogramPdf::normalize() -> fptype {
    Smoother smoother;
    smoother.funcIdx = getFunctionIndex();

    thrust::counting_iterator<int> binIndex(0);
    thrust::transform(binIndex, binIndex + dev_base_histogram->size(), dev_smoothed_histogram->begin(), smoother);

    // return totalEvents;
    fptype ret = thrust::reduce(dev_smoothed_histogram->begin(), dev_smoothed_histogram->end());

    for(unsigned int varIndex = 0; varIndex < observablesList.size(); ++varIndex) {
        fptype binSize = constantsList[3 + 3 * varIndex + 2];
        ret *= binSize; // Bin size cached by constructor.
    }

    // if (cpuDebug & 1) std::cout << "Normalizing " << getName() << " " << host_params[host_indices[parameters + 1]] <<
    // " " << ret << std::endl;
    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    return ret;
}
} // namespace GooFit
