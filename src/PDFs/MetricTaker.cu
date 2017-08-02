#include "goofit/PDFs/MetricTaker.h"

#include "goofit/GlobalCudaDefines.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/detail/ThrustOverride.h"

#include "goofit/BinnedDataSet.h"
#include "goofit/Error.h"
#include "goofit/FitControl.h"
#include "goofit/Log.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/Variable.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>

namespace GooFit {

__device__ fptype MetricTaker::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculate event offset for this thread.
    int eventIndex       = thrust::get<0>(t);
    int eventSize        = thrust::get<2>(t);
    fptype *eventAddress = thrust::get<1>(t) + (eventIndex * abs(eventSize));

    // Causes stack size to be statically undeterminable.
    fptype ret = callFunction(eventAddress, functionIdx, parameters);

    // Notice assumption here! For unbinned fits the 'eventAddress' pointer won't be used
    // in the metric, so it doesn't matter what it is. For binned fits it is assumed that
    // the structure of the event is (obs1 obs2... binentry binvolume), so that the array
    // passed to the metric consists of (binentry binvolume).
    ret = (*(reinterpret_cast<device_metric_ptr>(device_function_table[metricIndex])))(
        ret, eventAddress + (abs(eventSize) - 2), parameters);
    return ret;
}

#define MAX_NUM_OBSERVABLES 5

__device__ fptype binCenters[1024*MAX_NUM_OBSERVABLES];

__device__ fptype MetricTaker::operator()(thrust::tuple<int, int, fptype*> t) const {
    // Bin index, event size, base address [lower, upper,getNumBins]

    int evtSize   = thrust::get<1>(t);
    int binNumber = thrust::get<0>(t);

    // To convert global bin number to (x,y,z...) coordinates: For each dimension, take the mod
    // with the number of bins in that dimension. Then divide by the number of bins, in effect
    // collapsing so the grid has one fewer dimension. Rinse and repeat.
    unsigned int *indices = paramIndices + parameters;

    for(int i = 0; i < evtSize; ++i) {
        fptype lowerBound = thrust::get<2>(t)[3 * i + 0];
        fptype upperBound = thrust::get<2>(t)[3 * i + 1];
        auto numBins      = static_cast<int>(floor(thrust::get<2>(t)[3 * i + 2] + 0.5));
        int localBin      = binNumber % numBins;

        fptype x = upperBound - lowerBound;
        x /= numBins;
        x *= (localBin + 0.5);
        x += lowerBound;
        binCenters[indices[indices[0] + 2 + i] + THREADIDX * MAX_NUM_OBSERVABLES] = x;
        binNumber /= numBins;
    }

    // Causes stack size to be statically undeterminable.
    fptype ret = callFunction(binCenters + THREADIDX * MAX_NUM_OBSERVABLES, functionIdx, parameters);
    return ret;
}

MetricTaker::MetricTaker(PdfBase *dat, void *dev_functionPtr)
    : metricIndex(0)
    , functionIdx(dat->getFunctionIndex())
    , parameters(dat->getParameterIndex()) {
    // std::cout << "MetricTaker constructor with " << functionIdx << std::endl;

    auto localPos = functionAddressToDeviceIndexMap.find(dev_functionPtr);

    if(localPos != functionAddressToDeviceIndexMap.end()) {
        metricIndex = (*localPos).second;
    } else {
        metricIndex                                      = num_device_functions;
        host_function_table[num_device_functions]        = dev_functionPtr;
        functionAddressToDeviceIndexMap[dev_functionPtr] = num_device_functions;
        num_device_functions++;
        MEMCPY_TO_SYMBOL(device_function_table,
                         host_function_table,
                         num_device_functions * sizeof(void *),
                         0,
                         cudaMemcpyHostToDevice);
    }
}

MetricTaker::MetricTaker(int fIdx, int pIdx)
    : metricIndex(0)
    , functionIdx(fIdx)
    , parameters(pIdx) {
    // This constructor should only be used for binned evaluation, ie for integrals.
}

} // namespace GooFit
