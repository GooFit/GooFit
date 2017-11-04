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

__device__ fptype device_EvalHistogram(fptype *evt, ParameterContainer &pc) {
    // Structure is
    // nP smoothingIndex totalHistograms (limit1 step1 bins1) (limit2 step2 bins2) nO o1 o2
    // where limit and step are indices into functorConstants.

    int numCons          = RO_CACHE(pc.constants[pc.constantIdx]);
    int numObs           = RO_CACHE(pc.observables[pc.observableIdx]);
    int numParms         = RO_CACHE(pc.parameters[pc.parameterIdx]);
    int numVars          = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    int globalBinNumber  = 0;
    int previous         = 1;
    int myHistogramIndex = RO_CACHE(pc.constants[pc.constantIdx + 2]); // 1 only used for smoothing

    for(int i = 0; i < numVars; ++i) {
        int varIndex      = RO_CACHE(pc.observables[pc.observableIdx + 1 + i]);
        int lowerBoundIdx = 3 * (i + 1);
        // if (gpuDebug & 1) printf("[%i, %i] Smoothed: %i %i %i\n", BLOCKIDX, THREADIDX, i, varIndex,
        // indices[varIndex]);
        fptype currVariable = evt[varIndex];
        fptype lowerBound   = RO_CACHE(pc.constants[pc.constantIdx + i * 3 + 5]);
        fptype step         = RO_CACHE(pc.constants[pc.constantIdx + i * 3 + 6]);

        currVariable -= lowerBound;
        currVariable /= step;
        // if (gpuDebug & 1) printf("[%i, %i] Smoothed: %i %i %f %f %f %f\n", BLOCKIDX, THREADIDX, i, varIndex,
        // currVariable, lowerBound, step, evt[varIndex]);

        auto localBinNumber = static_cast<int>(floor(currVariable));
        globalBinNumber += previous * localBinNumber;

        // This is accessing too far ahead?
        int offset = RO_CACHE(pc.constants[pc.constantIdx + lowerBoundIdx + 1]);
        previous *= offset;
    }

    fptype *myHistogram = dev_smoothed_histograms[myHistogramIndex];
    fptype ret          = myHistogram[globalBinNumber];

    pc.incrementIndex(1, numParms, numCons, numObs, 1);

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
    int funcIdx;

    __device__ fptype operator()(int globalBin) {
        // unsigned int *indices = paramIndices + parameters;
        ParameterContainer pc;

        while(pc.funcIdx < funcIdx)
            pc.incrementIndex();

        int numObs       = RO_CACHE(pc.observables[pc.observableIdx]);
        int numVars      = RO_CACHE(pc.constants[pc.constantIdx + 1]);
        fptype smoothing = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
        // int histIndex         = RO_CACHE(indices[2]);
        // brad: need to debug this variable.
        int histIndex       = RO_CACHE(pc.constants[pc.constantIdx + 2]);
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
                // int lowerBoundIdx   = 3*(i+1);
                // int localNumBins = indices[6 + v*4];
                int localNumBins = RO_CACHE(pc.constants[pc.constantIdx + 3 * (v + 1) + 1]);
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
    int numVars = hist->numVariables();
    // int numConstants = 2 * numVars;
    // registerConstants(numConstants);
    // host_constants = new fptype[numConstants];
    totalEvents = 0;

    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(smoothing));
    pindices.push_back(totalHistograms);

    int varIndex = 0;

    constantsList.push_back(hist->getVariables().size());
    constantsList.push_back(totalHistograms);
    constantsList.push_back(numVars);

    for(Variable *var : hist->getVariables()) {
        registerObservable(var);
        // pindices.push_back((*var)->index);
        // pindices.push_back(cIndex + 2 * varIndex + 0);
        // pindices.push_back(cIndex + 2 * varIndex + 1);
        // pindices.push_back(var->getNumBins());
        constantsList.push_back(var->getNumBins());
        // constantsList.push_back (varIndex);
        constantsList.push_back(var->getLowerLimit());
        constantsList.push_back(var->getBinSize());

        // host_constants[2 * varIndex + 0] = var->getLowerLimit(); // NB, do not put cIndex here, it is accounted for
        // by
        // the offset in MEMCPY_TO_SYMBOL below.
        // host_constants[2 * varIndex + 1] = var->getBinSize();
        varIndex++;
    }

    unsigned int numbins = hist->getNumBins();
    thrust::host_vector<fptype> host_histogram;

    for(unsigned int i = 0; i < numbins; ++i) {
        fptype curr = hist->getBinContent(i);
        host_histogram.push_back(curr);
        totalEvents += curr;
    }

    // MEMCPY_TO_SYMBOL(functorConstants,
    //                 host_constants,
    //                 numConstants * sizeof(fptype),
    //                 cIndex * sizeof(fptype),
    //                 cudaMemcpyHostToDevice);

    if(totalEvents > 0)
        copyHistogramToDevice(host_histogram);
    else
        std::cout << "Warning: Empty histogram supplied to " << getName()
                  << " not copied to device. Expect copyHistogramToDevice call later.\n";

    GET_FUNCTION_ADDR(ptr_to_EvalHistogram);
    initialize(pindices);
}

void SmoothHistogramPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_EvalHistogram);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_EvalHistogram");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
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

    for(auto &observable : observablesList) {
        expectedBins *= observable->getNumBins();
    }

    if(expectedBins != host_histogram.size()) {
        std::cout << "Warning: Histogram supplied to " << getName() << " has " << host_histogram.size()
                  << " bins, expected " << expectedBins << " - may indicate a problem.\n";
    }
}

__host__ fptype SmoothHistogramPdf::normalize() const {
    Smoother smoother;
    smoother.funcIdx = getFunctionIndex();

    thrust::counting_iterator<int> binIndex(0);
    thrust::transform(binIndex, binIndex + dev_base_histogram->size(), dev_smoothed_histogram->begin(), smoother);

    // return totalEvents;
    fptype ret = thrust::reduce(dev_smoothed_histogram->begin(), dev_smoothed_histogram->end());

    int obs = constantsList[0];
    for(unsigned int varIndex = 0; varIndex < observablesList.size(); ++varIndex) {
        fptype binSize = constantsList[obs + 3 * varIndex + 3];
        ret *= binSize; // Bin size cached by constructor.
    }

    // if (cpuDebug & 1) std::cout << "Normalising " << getName() << " " << host_params[host_indices[parameters + 1]] <<
    // " " << ret << std::endl;
    host_normalisations[normalIdx + 1] = 1.0 / ret;
    return ret;
}
} // namespace GooFit
