#include <goofit/PDFs/MetricTaker.h>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/detail/ThrustOverride.h>

#include <goofit/BinnedDataSet.h>
#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/Log.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

namespace GooFit {

__device__ fptype MetricTaker::operator()(thrust::tuple<int, fptype *, int> t) const {
    ParameterContainer pc;

    // Calculate event offset for this thread.
    int eventIndex       = thrust::get<0>(t);
    int eventSize        = thrust::get<2>(t);
    fptype *eventAddress = thrust::get<1>(t) + (eventIndex * abs(eventSize));

    int idx = abs(eventSize - 2);
    if(idx < 0)
        idx = 0;

    // fptype obs  = events[idx];
    fptype norm = pc.getNormalization(0);

    // Causes stack size to be statically undeterminable.
    fptype ret = callFunction(eventAddress, pc);

    // Notice assumption here! For unbinned fits the 'eventAddress' pointer won't be used
    // in the metric, so it doesn't matter what it is. For binned fits it is assumed that
    // the structure of the event is (obs1 obs2... binentry binvolume), so that the array
    // passed to the metric consists of (binentry binvolume).

    ret = (*(reinterpret_cast<device_metric_ptr>(d_function_table[pc.funcIdx])))(
        ret, eventAddress + (abs(eventSize) - 2), norm);

    return ret;
}

#define MAX_NUM_OBSERVABLES 10

__device__ fptype MetricTaker::operator()(thrust::tuple<int, int, fptype *> t) const {
    ParameterContainer pc;

    int evtSize   = thrust::get<1>(t);
    int binNumber = thrust::get<0>(t);

    auto *events = new fptype[10];

    // for (int i = 0; i < evtSize; i++)
    //    pc.events[i] =

    // Do not understand why this cannot be declared __shared__. Dynamically allocating shared memory is apparently
    // complicated.
    // fptype* binCenters = (fptype*) malloc(evtSize * sizeof(fptype));
    // fptype binCenters[MAX_NUM_OBSERVABLES];

    // To convert global bin number to (x,y,z...) coordinates: For each dimension, take the mod
    // with the number of bins in that dimension. Then divide by the number of bins, in effect
    // collapsing so the grid has one fewer dimension. Rinse and repeat.

    // before we start, we need to progress ourselves based on the function we need to run
    for(int i = 0; i < functionIdx; i++)
        pc.incrementIndex(); // need to use the slow version, since we are not starting from index 0.

    // put our index here...

    for(int i = 0; i < evtSize; ++i) {
        int id = pc.getObservable(i);
        fptype lowerBound = thrust::get<2>(t)[3 * i + 0];
        fptype upperBound = thrust::get<2>(t)[3 * i + 1];
        auto numBins      = static_cast<int>(floor(thrust::get<2>(t)[3 * i + 2] + 0.5));
        int localBin      = binNumber % numBins;

        fptype x = upperBound - lowerBound;
        x /= numBins;
        x *= (localBin + 0.5);
        x += lowerBound;
        // binCenters[i + THREADIDX*MAX_NUM_OBSERVABLES] = x;
        events[id] = x;
        binNumber /= numBins;
    }

    // Causes stack size to be statically undeterminable.
    fptype ret = callFunction(events, pc);

    delete[] events;

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
        metricIndex                                      = host_function_table.size();
        functionAddressToDeviceIndexMap[dev_functionPtr] = host_function_table.size();
        host_function_table.push_back(dev_functionPtr);
        host_function_table.sync(d_function_table);
    }
}

MetricTaker::MetricTaker(int fIdx, int pIdx)
    : metricIndex(0)
    , functionIdx(fIdx)
    , parameters(pIdx) {
    // This constructor should only be used for binned evaluation, ie for integrals.
}

} // namespace GooFit
