#include "goofit/PDFs/basic/SmoothHistogramPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

__constant__ fptype *dev_base_histograms[100]; // Multiple histograms for the case of multiple PDFs
__constant__ fptype *dev_smoothed_histograms[100];
unsigned int SmoothHistogramPdf::totalHistograms = 0;

__device__ int dev_powi(int base, int exp) {
    int ret = 1;

    for(int i = 0; i < exp; ++i)
        ret *= base;

    return ret;
}

__device__ fptype device_EvalHistogram(fptype *evt, fptype *p, unsigned int *indices) {
    // Structure is
    // nP smoothingIndex totalHistograms (limit1 step1 bins1) (limit2 step2 bins2) nO o1 o2
    // where limit and step are indices into functorConstants.

    int numVars          = RO_CACHE(indices[RO_CACHE(indices[0]) + 1]);
    int globalBinNumber  = 0;
    int previous         = 1;
    int myHistogramIndex = RO_CACHE(indices[2]); // 1 only used for smoothing

    for(int i = 0; i < numVars; ++i) {
        int varIndex      = RO_CACHE(indices[RO_CACHE(indices[0]) + 2 + i]);
        int lowerBoundIdx = 3 * (i + 1);
        // if (gpuDebug & 1) printf("[%i, %i] Smoothed: %i %i %i\n", BLOCKIDX, THREADIDX, i, varIndex,
        // indices[varIndex]);
        fptype currVariable = evt[varIndex];
        fptype lowerBound   = RO_CACHE(functorConstants[RO_CACHE(indices[lowerBoundIdx + 0])]);
        fptype step         = RO_CACHE(functorConstants[RO_CACHE(indices[lowerBoundIdx + 1])]);

        currVariable -= lowerBound;
        currVariable /= step;
        // if (gpuDebug & 1) printf("[%i, %i] Smoothed: %i %i %f %f %f %f\n", BLOCKIDX, THREADIDX, i, varIndex,
        // currVariable, lowerBound, step, evt[varIndex]);

        auto localBinNumber = static_cast<int>(floor(currVariable));
        globalBinNumber += previous * localBinNumber;
        previous *= indices[lowerBoundIdx + 2];
    }

    fptype *myHistogram = dev_smoothed_histograms[myHistogramIndex];
    fptype ret          = RO_CACHE(myHistogram[globalBinNumber]);

    // if ((gpuDebug & 1) && (evt[8] < 0.5) && (paramIndices + debugParamIndex == indices)) printf("Smoothed: %f %f %f
    // %i %f\n", evt[6], evt[7], myHistogram[globalBinNumber], globalBinNumber,
    // dev_base_histograms[myHistogramIndex][globalBinNumber]);
    // if (gpuDebug & 1) printf("Smoothed: %f %f %f %i %f\n", evt[0], evt[1], myHistogram[globalBinNumber],
    // globalBinNumber, dev_base_histograms[myHistogramIndex][globalBinNumber]);
    // if (gpuDebug & 1) printf("Smoothed: %f %f %f %i %f %f\n", evt[0], evt[1], ret, globalBinNumber,
    // dev_base_histograms[myHistogramIndex][globalBinNumber], p[indices[1]]);
    return ret;
}

struct Smoother {
    int parameters;

    __device__ fptype operator()(int globalBin) {
        unsigned int *indices = paramIndices + parameters;
        int numVars           = RO_CACHE(indices[RO_CACHE(indices[0]) + 1]);
        fptype smoothing      = RO_CACHE(cudaArray[RO_CACHE(indices[1])]);
        int histIndex         = RO_CACHE(indices[2]);
        fptype *myHistogram   = dev_base_histograms[histIndex];
        fptype centralValue   = myHistogram[globalBin];

        fptype otherBinsTotal  = 0;
        int numSurroundingBins = 0;
        int otherBins          = dev_powi(3, numVars);

        for(int i = 0; i < otherBins; ++i) {
            int currBin       = globalBin;
            int localPrevious = 1;
            int trackingBin   = globalBin;
            bool offSomeAxis  = false;

            for(int v = 0; v < numVars; ++v) {
                // int lowerBoundIdx   = 3*(i+1);
                // int localNumBins = indices[6 + v*4];
                int localNumBins = RO_CACHE(indices[3 * (v + 1) + 2]);
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

        // if (7010 == globalBin) printf("Smoothing: %f %f %f %i %f\n", myHistogram[globalBin], otherBinsTotal,
        // smoothing, numSurroundingBins, centralValue);
        return centralValue;
    }
};

__device__ device_function_ptr ptr_to_EvalHistogram = device_EvalHistogram;

__host__ SmoothHistogramPdf::SmoothHistogramPdf(std::string n, BinnedDataSet *hist, Variable *smoothing)
    : GooPdf(nullptr, n) {
    int numVars      = hist->numVariables();
    int numConstants = 2 * numVars;
    registerConstants(numConstants);
    host_constants = new fptype[numConstants];
    totalEvents    = 0;

    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(smoothing));
    pindices.push_back(totalHistograms);

    int varIndex = 0;

    for(Variable *var : hist->getVariables()) {
        registerObservable(var);
        // pindices.push_back((*var)->index);
        pindices.push_back(cIndex + 2 * varIndex + 0);
        pindices.push_back(cIndex + 2 * varIndex + 1);
        pindices.push_back(var->getNumBins());

        host_constants[2 * varIndex + 0] = var->getLowerLimit(); // NB, do not put cIndex here, it is accounted for by
                                                                 // the offset in MEMCPY_TO_SYMBOL below.
        host_constants[2 * varIndex + 1] = var->getBinSize();
        varIndex++;
    }

    unsigned int numbins = hist->getNumBins();
    thrust::host_vector<fptype> host_histogram;

    for(unsigned int i = 0; i < numbins; ++i) {
        fptype curr = hist->getBinContent(i);
        host_histogram.push_back(curr);
        totalEvents += curr;
    }

    MEMCPY_TO_SYMBOL(functorConstants,
                     host_constants,
                     numConstants * sizeof(fptype),
                     cIndex * sizeof(fptype),
                     cudaMemcpyHostToDevice);

    if(totalEvents > 0)
        copyHistogramToDevice(host_histogram);
    else
        std::cout << "Warning: Empty histogram supplied to " << getName()
                  << " not copied to device. Expect copyHistogramToDevice call later.\n";

    GET_FUNCTION_ADDR(ptr_to_EvalHistogram);
    initialize(pindices);
}

fptype *pointerToFirst(thrust::device_vector<fptype> *hist) { return (&((*hist)[0])).get(); }

fptype *pointerToFirst(thrust::host_vector<fptype> *hist) {
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

    for(auto &observable : observables) {
        expectedBins *= observable->getNumBins();
    }

    if(expectedBins != host_histogram.size()) {
        std::cout << "Warning: Histogram supplied to " << getName() << " has " << host_histogram.size()
                  << " bins, expected " << expectedBins << " - may indicate a problem.\n";
    }
}

__host__ fptype SmoothHistogramPdf::normalize() const {
    Smoother smoother;
    smoother.parameters = parameters;

    thrust::counting_iterator<int> binIndex(0);
    thrust::transform(binIndex, binIndex + dev_base_histogram->size(), dev_smoothed_histogram->begin(), smoother);

    // return totalEvents;
    fptype ret = thrust::reduce(dev_smoothed_histogram->begin(), dev_smoothed_histogram->end());

    for(unsigned int varIndex = 0; varIndex < observables.size(); ++varIndex) {
        ret *= host_constants[2 * varIndex + 1]; // Bin size cached by constructor.
    }

    // if (cpuDebug & 1) std::cout << "Normalising " << getName() << " " << host_params[host_indices[parameters + 1]] <<
    // " " << ret << std::endl;
    host_normalisation[parameters] = 1.0 / ret;
    return ret;
}
} // namespace GooFit
