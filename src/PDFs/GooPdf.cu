#include <goofit/BinnedDataSet.h>
#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/Log.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/Version.h>
#include <goofit/detail/SmartVector.h>
#include <goofit/detail/ThrustOverride.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

// These variables are either function-pointer related (thus specific to this implementation)
// or constrained to be in the CUDAglob translation unit by nvcc limitations; otherwise they
// would be in PdfBase.

// For debugging

__constant__ int gpuDebug;
__constant__ unsigned int debugParamIndex;

int cpuDebug = 0;

GooPdf::~GooPdf() { cleanup(); }

// Reduce the PDFs to a single value based on metric taker
// numVars will be different for binned or unbinned fit
// This does NOT normalize!
__host__ double GooPdf::reduce_with_metric() const {
    double ret;

    double start = 0.0;

    thrust::constant_iterator<int> eventSize(get_event_size());
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);

#ifdef GOOFIT_MPI
    size_t entries_to_process = m_iEventsPerTask;
#else
    size_t entries_to_process = numEntries;
#endif

    // Calls and sums in parallel:
    // logger(0, arrayAddress, eventSize) +
    // logger(1, arrayAddress, eventSize) +
    // ...

    ret = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + entries_to_process, arrayAddress, eventSize)),
        *logger,
        start,
        thrust::plus<double>());

#ifdef GOOFIT_MPI
    double r = ret;
    MPI_Allreduce(&r, &ret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    return ret;
}

// Reduce the PDFs to a single value based on metric taker
// numVars will be different for binned or unbinned fit
// This does NOT normalize!
__host__ double GooPdf::reduce_with_bins() const {
    double sum;
    double start = 0.0;

    thrust::constant_iterator<int> eventSize(observablesList.size());
    thrust::constant_iterator<fptype *> arrayAddress(normRanges);
    thrust::counting_iterator<int> binIndex(0);
    logger->setFunctionIndex(functionIdx);

    size_t totalBins = get_bin_grid_size();

    // Calls and sums in parallel:
    // logger(0, eventSize, arrayAddress) +
    // logger(1, eventSize, arrayAddress) +
    // ...
    // Note that these are *reversed*!

    sum = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
        thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
        *logger,
        start,
        thrust::plus<double>());

    return sum;
}

/// This evaluates the current function over the data. Does *not* prepare
/// or normalize
__host__ void GooPdf::evaluate_with_metric(thrust::device_vector<fptype> &results) const {
    // if(results.size() != numEntries)
    //    results.resize(numEntries);

    thrust::constant_iterator<int> eventSize(observablesList.size());
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);

#ifdef GOOFIT_MPI
    size_t entries_to_process = m_iEventsPerTask;
#else
    size_t entries_to_process = numEntries;
#endif

    results.resize(entries_to_process);

    // Calls in parallel:
    // logger(0, arrayAddress, eventSize)
    // logger(1, arrayAddress, eventSize)
    // ...
    thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + entries_to_process, arrayAddress, eventSize)),
        results.begin(),
        *logger);

    // We need to copy each 'results' buffer to each other
#ifdef GOOFIT_MPI
    // copy our local device buffer to a temporary host_vector
    thrust::host_vector<fptype> local_results = results;

    int myId, numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    int counts[numProcs];
    int displs[numProcs];

    // gather all our counts.
    MPI_Allgather(&entries_to_process, 1, MPI_INT, &counts[0], 1, MPI_INT, MPI_COMM_WORLD);

    // calculate our displs.
    displs[0] = 0;
    for(int i = 1; i < numProcs; i++)
        displs[i] = displs[i - 1] + counts[i - 1];

    thrust::host_vector<fptype> total_results;
    total_results.resize(numEntries);

    MPI_Allgatherv(&local_results[0],
                   local_results.size(),
                   MPI_DOUBLE,
                   &total_results[0],
                   &counts[0],
                   &displs[0],
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);

    // copy our results back to our device_vector.
    results = total_results;
#endif
}

__host__ thrust::host_vector<fptype> GooPdf::evaluate_with_metric() const {
    thrust::device_vector<fptype> results;
    results.resize(numEntries);
    evaluate_with_metric(results);
    return thrust::host_vector<fptype>(results);
}

__host__ void GooPdf::setIndices() {
    // If not set, perform unbinned Nll fit!
    if(!fitControl) {
        setFitControl(std::make_shared<UnbinnedNllFit>());
    }

    // Ensure that we properly populate *logger with the correct metric
    setMetrics();

    PdfBase::setIndices();

    GOOFIT_DEBUG("host_function_table[{}] = {} (fitControl)", host_function_table.size(), fitControl->getName());

    host_function_table.push_back(getMetricPointer(fitControl->getMetric()));
    host_function_name.push_back(fitControl->getName());

    // copy all the device functions over:
    pre_run();
}

__host__ int GooPdf::findFunctionIdx(void *dev_functionPtr) {
    // Code specific to function-pointer implementation
    auto localPos = functionAddressToDeviceIndexMap.find(dev_functionPtr);

    if(localPos != functionAddressToDeviceIndexMap.end()) {
        return (*localPos).second;
    }

    int fIdx = host_function_table.size();
    host_function_table.push_back(dev_functionPtr);
    host_function_name.push_back("findFunctionIdx");
    functionAddressToDeviceIndexMap[dev_functionPtr] = fIdx;
    host_function_table.sync(d_function_table);

    return fIdx;
}

/*__host__ auto GooPdf::lookUpFunctionIdx(void *dev_functionPtr) -> int {
    // 211218 mds
    // derived from findFunctionIdx, but only looks up value
    // returns -99 if not found; does not add anything to
    // host_function_table or host_function_name

    // Code specific to function-pointer implementation
    auto localPos = functionAddressToDeviceIndexMap.find(dev_functionPtr);

    if(localPos != functionAddressToDeviceIndexMap.end()) {
        return (*localPos).second;
    }

    int fIdx = -99;

    return fIdx;
}
*/
__host__ void GooPdf::setDebugMask(int mask, bool setSpecific) const {
    cpuDebug = mask;
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
    gpuDebug = cpuDebug;

    if(setSpecific)
        debugParamIndex = parameters;

#else
    MEMCPY_TO_SYMBOL(gpuDebug, &cpuDebug, sizeof(int), 0, cudaMemcpyHostToDevice);

    if(setSpecific)
        MEMCPY_TO_SYMBOL(debugParamIndex, &parameters, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);

#endif
}

__host__ double GooPdf::calculateNLL() {
    GOOFIT_MAYBE_UNUSED fptype norm = normalize();
    GOOFIT_TRACE("GooPdf::calculateNLL calling normalize: {} (host_norm should be 1: {})",
                 norm,
                 host_normalizations.at(normalIdx + 1));

    if(host_normalizations[normalIdx + 1] <= 0)
        GooFit::abort(__FILE__, __LINE__, getName() + " non-positive normalization", this);

    // make this memcpy async
    host_normalizations.sync(d_normalizations);

    fptype ret = reduce_with_metric();
    GOOFIT_TRACE("GooPdf::calculateNLL calling sumOfNll = {}", ret);

    if(0.0 == ret)
        GooFit::abort(__FILE__, __LINE__, getName() + " zero NLL", this);

    return 2.0 * ret;
}

__host__ std::vector<fptype> GooPdf::evaluateAtPoints(Observable var) {
    setFitControl(std::make_shared<EvalFit>());

    setIndices();

    normalize();
    host_normalizations.sync(d_normalizations);

    UnbinnedDataSet tempdata(observablesList);

    double step = var.getBinSize();

    for(int i = 0; i < var.getNumBins(); ++i) {
        var.setValue(var.getLowerLimit() + (i + 0.5) * step);
        tempdata.addEvent();
    }

    auto old = getData();
    setData(&tempdata);

    normalize();

    host_normalizations.sync(d_normalizations);

    thrust::host_vector<fptype> h_results = evaluate_with_metric();
    std::vector<fptype> res;
    res.resize(var.getNumBins());

    for(int i = 0; i < var.getNumBins(); ++i) {
        fptype n = host_normalizations.at(normalIdx + 1);
        fptype v = h_results[i];
        res[i]   = v * n;
    }

    setData(old);

    return res;
}

__host__ fptype GooPdf::getValue(EvalFunc evalfunc) {
    if(evalfunc == EvalFunc::Prob) {
        setFitControl(std::make_shared<ProbFit>());
    } else if(evalfunc == EvalFunc::Eval) {
        setFitControl(std::make_shared<EvalFit>());
    } else if(evalfunc == EvalFunc::NLL) {
        setFitControl(std::make_shared<UnbinnedNllFit>());
    } else
        throw GeneralError("That EvalFunc is not supported");

    setIndices();
    normalize();
    host_normalizations.sync(d_normalizations);

    UnbinnedDataSet point(observablesList);
    point.addEvent();
    auto old = getData();
    setData(&point);

    auto results = evaluate_with_metric();

    setData(old);

    return results[0];
}

__host__ fptype GooPdf::normalize() {
    if(!fitControl->metricIsPdf()) {
        GOOFIT_TRACE("{}: metricIsPdf, returning 1", getName());
        host_normalizations.at(normalIdx + 1) = 1.0;
        cachedNormalization                   = 1.0;
        return 1.0;
    }

    fptype ret = 1;

    if(hasAnalyticIntegral()) {
        // Loop goes only over observables of this PDF.
        for(const Observable &v : observablesList) {
            GOOFIT_TRACE("{}: Analytically integrating over {}", getName(), v.getName());
            ret *= integrate(v.getLowerLimit(), v.getUpperLimit());
        }

        host_normalizations.at(normalIdx + 1) = 1.0 / ret;
        cachedNormalization                   = 1.0 / ret;
        GOOFIT_TRACE("{}: Param {} integral is = {}", getName(), parameters, ret);

        return ret;
    }

    GOOFIT_TRACE("{}, Computing integral without analytic help", getName());

    size_t totalBins = get_bin_grid_size();

    GOOFIT_TRACE("Total Bins: {}", totalBins);

    ret *= get_bin_grid_volume();
    GOOFIT_TRACE("Bin volume: {}", ret);

    ret /= totalBins;

    GOOFIT_TRACE("Total bin volume: {}", ret);

    fptype sum = reduce_with_bins();

    GOOFIT_TRACE("Sum: {}", sum);

    if(std::isnan(sum)) {
        GooFit::abort(__FILE__, __LINE__, getName() + " NaN in normalization", this);
    } else if(0 >= sum) {
        GooFit::abort(__FILE__, __LINE__, "Non-positive normalization (mds) ", this);
    }

    ret *= sum;

    if(0 == ret)
        GooFit::abort(__FILE__, __LINE__, "Zero integral");

    GOOFIT_TRACE("{}: Param {} integral is ~= {}", getName(), normalIdx, ret);
    host_normalizations.at(normalIdx + 1) = 1.0 / ret;
    cachedNormalization                   = 1.0 / ret;
    return (fptype)ret;
}

__device__ double callFunction(fptype *eventAddress, ParameterContainer &pc) {
    return (*(reinterpret_cast<device_function_ptr>(d_function_table[pc.funcIdx])))(eventAddress, pc);
}

__host__ std::vector<std::vector<fptype>> GooPdf::getCompProbsAtDataPoints() {
    // note, we need to overwrite what our metric operator is going to do, and restore previous
    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());
    normalize();
    host_normalizations.sync(d_normalizations);
    thrust::host_vector<fptype> host_results = evaluate_with_metric();

    std::vector<std::vector<fptype>> values;
    values.resize(components.size() + 1);
    values[0] = std::vector<fptype>(host_results.begin(), host_results.end());

    for(unsigned int i = 0; i < components.size(); ++i) {
        // we need to recreate the indexing for each component
        // components[i]->setFitControl(std::make_shared<ProbFit>());
        components[i]->setIndices();
        components[i]->normalize();

        GOOFIT_TRACE("host_function_table[{}] = {}", host_function_table.size(), fitControl->getName());
        host_function_table.push_back(getMetricPointer(fitControl->getMetric()));
        host_function_name.push_back(fitControl->getName());

        // copy all the device functions over:
        GOOFIT_DEBUG("Copying all host side parameters to device (normalizations too)");
        pre_run();

        auto result   = evaluate_with_metric();
        values[1 + i] = std::vector<fptype>(result.begin(), result.end());
    }

    // restore previous fit control
    setFitControl(fc);
    return values;
}
/*
__host__ void GooPdf::listAllComponents(PdfBase* someComponent) {

   auto moreComponents = someComponent->getComponents();
   for (auto &mComponent : moreComponents) {
     std::cout << " mComponent = " << mComponent << std::endl;
     std::cout << "  *mComponent = " << *mComponent << std::endl;
     listAllComponents(mComponent);
   }
}
*/
} // namespace GooFit
