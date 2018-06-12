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

__host__ void GooPdf::setIndices() {
    // If not set, perform unbinned Nll fit!
    if(!fitControl)
        setFitControl(std::make_shared<UnbinnedNllFit>());

    // Ensure that we properly populate *logger with the correct metric
    setMetrics();

    GOOFIT_TRACE("GooPdf::setIndices!");
    PdfBase::setIndices();

    GOOFIT_DEBUG("host_function_table[{}] = {} (fitControl)", num_device_functions, fitControl->getName());
    if(num_device_functions >= GOOFIT_MAXFUNC)
        throw GeneralError("Too many device functions! Set GOOFIT_MAXFUNC to a larger value than {}", GOOFIT_MAXFUNC);
    host_function_table[num_device_functions++] = getMetricPointer(fitControl->getMetric());

    // copy all the device functions over:
    GOOFIT_TRACE("Copying all host side parameters to device");
    MEMCPY_TO_SYMBOL(
        device_function_table, &host_function_table, num_device_functions * sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(d_parameters, &host_parameters, totalParameters * sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(d_constants, &host_constants, totalConstants * sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(d_observables, &host_observables, totalObservables * sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(
        d_normalizations, &host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);
}

__host__ int GooPdf::findFunctionIdx(void *dev_functionPtr) {
    // Code specific to function-pointer implementation
    auto localPos = functionAddressToDeviceIndexMap.find(dev_functionPtr);

    if(localPos != functionAddressToDeviceIndexMap.end()) {
        return (*localPos).second;
    }

    int fIdx = num_device_functions;
    if(num_device_functions >= GOOFIT_MAXFUNC)
        throw GeneralError("Too many device functions! Set GOOFIT_MAXFUNC to a larger value than {}", GOOFIT_MAXFUNC);
    host_function_table[num_device_functions]        = dev_functionPtr;
    functionAddressToDeviceIndexMap[dev_functionPtr] = num_device_functions;
    num_device_functions++;
    MEMCPY_TO_SYMBOL(
        device_function_table, host_function_table, num_device_functions * sizeof(void *), 0, cudaMemcpyHostToDevice);

#ifdef PROFILING
    host_timeHist[fIdx] = 0;
    MEMCPY_TO_SYMBOL(timeHistogram, host_timeHist, 10000 * sizeof(fptype), 0);
#endif

    return fIdx;
}

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

__host__ double GooPdf::sumOfNll(int numVars) const {
    static thrust::plus<double> cudaPlus;
    thrust::constant_iterator<int> eventSize(numVars);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);
    double dummy = 0;

    double ret;
#ifdef GOOFIT_MPI
    double r = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);

    MPI_Allreduce(&r, &ret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    ret = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#endif
    return ret;
}

__host__ double GooPdf::calculateNLL() {
    GOOFIT_MAYBE_UNUSED fptype norm = normalize();
    GOOFIT_TRACE("GooPdf::calculateNLL calling normalize: {} (host_norm should be 1: {})",
                 norm,
                 host_normalizations[normalIdx + 1]);

    if(host_normalizations[normalIdx + 1] <= 0)
        GooFit::abort(__FILE__, __LINE__, getName() + " non-positive normalization", this);

    // make this memcpy async
    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    auto numVars = (int)observablesList.size();

    if(fitControl->binnedFit()) {
        numVars += 2;
        numVars *= -1;
    }

    fptype ret = sumOfNll(numVars);
    GOOFIT_TRACE("GooPdf::calculateNLL calling sumOfNll = {}", ret);

    if(0.0 == ret)
        GooFit::abort(__FILE__, __LINE__, getName() + " zero NLL", this);

    return 2.0 * ret;
}

__host__ std::vector<fptype> GooPdf::evaluateAtPoints(Observable var) {
    setFitControl(std::make_shared<EvalFit>());

    setIndices();

    normalize();

    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);
    UnbinnedDataSet tempdata(observablesList);

    double step = var.getBinSize();

    for(int i = 0; i < var.getNumBins(); ++i) {
        var.setValue(var.getLowerLimit() + (i + 0.5) * step);
        tempdata.addEvent();
    }

    auto old = getData();
    setData(&tempdata);

    normalize();

    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    thrust::counting_iterator<int> eventIndex(0);
    thrust::constant_iterator<int> eventSize(observablesList.size());
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::device_vector<fptype> results(var.getNumBins());

#ifdef GOOFIT_MPI
    thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
        results.begin(),
        *logger);
#else
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                      results.begin(),
                      *logger);
#endif

    // Note, This is not fully realized with MPI.  We need to copy each 'results' buffer to each other 'MPI_Scatterv',
    // then we can do the rest.
    thrust::host_vector<fptype> h_results = results;
    std::vector<fptype> res;
    res.resize(var.getNumBins());

    for(int i = 0; i < var.getNumBins(); ++i) {
        fptype n = host_normalizations[normalIdx + 1];
        fptype v = h_results[i];
        res[i]   = v * n;
    }

    setData(old);

    return res;
}

__host__ fptype GooPdf::getValue(EvalFunc evalfunc) {
    if(evalfunc == EvalFunc::Prob)
        setFitControl(std::make_shared<ProbFit>());
    else if(evalfunc == EvalFunc::Eval)
        setFitControl(std::make_shared<EvalFit>());
    else if(evalfunc == EvalFunc::NLL)
        setFitControl(std::make_shared<UnbinnedNllFit>());
    else
        throw GeneralError("That EvalFunc is not supported");

    setIndices();
    normalize();
    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    UnbinnedDataSet point(observablesList);
    point.addEvent();
    auto old = getData();
    setData(&point);

    thrust::counting_iterator<int> eventIndex(0);
    thrust::constant_iterator<int> eventSize(observablesList.size());
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::device_vector<fptype> results(1);

    // MetricTaker evalor(this, getMetricPointer(evalfunc));
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + 1, arrayAddress, eventSize)),
                      results.begin(),
                      *logger);

    setData(old);

    return results[0];
}

__host__ fptype GooPdf::normalize() {
    if(!fitControl->metricIsPdf()) {
        GOOFIT_TRACE("{}: metricIsPdf, returning 1", getName());
        host_normalizations[normalIdx + 1] = 1.0;
        cachedNormalization                = 1.0;
        return 1.0;
    }

    fptype ret = 1;

    if(hasAnalyticIntegral()) {
        // Loop goes only over observables of this PDF.
        for(const Observable &v : observablesList) {
            GOOFIT_TRACE("{}: Analytically integrating over {}", getName(), v.getName());
            ret *= integrate(v.getLowerLimit(), v.getUpperLimit());
        }

        host_normalizations[normalIdx + 1] = 1.0 / ret;
        cachedNormalization                = 1.0 / ret;
        GOOFIT_TRACE("{}: Param {} integral is = {}", getName(), parameters, ret);

        return ret;
    }

    GOOFIT_TRACE("{}, Computing integral without analytic help", getName());

    int totalBins = 1;

    for(const Observable &v : observablesList) {
        ret *= v.getUpperLimit() - v.getLowerLimit();
        totalBins *= integrationBins > 0 ? integrationBins : v.getNumBins();

        GOOFIT_TRACE("Total bins {} due to {} {} {}", totalBins, v.getName(), integrationBins, v.getNumBins());
    }

    ret /= totalBins;

    fptype dummy = 0;
    static thrust::plus<fptype> cudaPlus;
    thrust::constant_iterator<fptype *> arrayAddress(normRanges);
    thrust::constant_iterator<int> eventSize(observablesList.size());
    thrust::counting_iterator<int> binIndex(0);
    // thrust::constant_iterator<int> funcIdx (functionIdx);
    logger->setFunctionIndex(functionIdx);

    fptype sum;
#ifdef GOOFIT_MPI
    fptype s = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
        thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
        *logger,
        dummy,
        cudaPlus);

    MPI_Allreduce(&s, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    sum = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
        thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
        *logger,
        dummy,
        cudaPlus);

    GOOFIT_TRACE("sum = {}", sum);
#endif

    if(std::isnan(sum)) {
        GooFit::abort(__FILE__, __LINE__, getName() + " NaN in normalization", this);
    } else if(0 >= sum) {
        GooFit::abort(__FILE__, __LINE__, "Non-positive normalization", this);
    }

    ret *= sum;

    if(0 == ret)
        GooFit::abort(__FILE__, __LINE__, "Zero integral");

    GOOFIT_TRACE("{}: Param {} integral is ~= {}", getName(), normalIdx, ret);
    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    return (fptype)ret;
}

__device__ fptype callFunction(fptype *eventAddress, ParameterContainer &pc) {
    return (*(reinterpret_cast<device_function_ptr>(device_function_table[pc.funcIdx])))(eventAddress, pc);
}

__host__ std::vector<std::vector<fptype>> GooPdf::getCompProbsAtDataPoints() {
    // note, we need to overwrite what our metric operator is going to do, and restore previous
    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());

    normalize();
    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    int numVars = observablesList.size();

    if(fitControl->binnedFit()) {
        numVars += 2;
        numVars *= -1;
    }

    thrust::device_vector<fptype> results(numEntries);
    thrust::constant_iterator<int> eventSize(numVars);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);
    // MetricTaker evalor(this, getMetricPointer("ptr_to_Prob"));
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                      results.begin(),
                      *logger);
    std::vector<std::vector<fptype>> values;
    values.resize(components.size() + 1);
    thrust::host_vector<fptype> host_results = results;

    for(unsigned int i = 0; i < host_results.size(); ++i) {
        values[0].push_back(host_results[i]);
    }

    for(unsigned int i = 0; i < components.size(); ++i) {
        // we need to recreate the indexing for each component
        // components[i]->setFitControl(std::make_shared<ProbFit>());
        components[i]->setIndices();
        components[i]->normalize();

        GOOFIT_TRACE("host_function_table[{}] = {}", num_device_functions, fitControl->getName());
        host_function_table[num_device_functions] = getMetricPointer(fitControl->getMetric());
        num_device_functions++;

        // copy all the device functions over:
        GOOFIT_DEBUG("Copying all host side parameters to device");
        MEMCPY_TO_SYMBOL(device_function_table,
                         &host_function_table,
                         num_device_functions * sizeof(fptype),
                         0,
                         cudaMemcpyHostToDevice);
        MEMCPY_TO_SYMBOL(d_parameters, &host_parameters, totalParameters * sizeof(fptype), 0, cudaMemcpyHostToDevice);
        MEMCPY_TO_SYMBOL(d_constants, &host_constants, totalConstants * sizeof(fptype), 0, cudaMemcpyHostToDevice);
        MEMCPY_TO_SYMBOL(
            d_observables, &host_observables, totalObservables * sizeof(fptype), 0, cudaMemcpyHostToDevice);
        MEMCPY_TO_SYMBOL(
            d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

        thrust::counting_iterator<int> ceventIndex(0);
        thrust::transform(
            thrust::make_zip_iterator(thrust::make_tuple(ceventIndex, arrayAddress, eventSize)),
            thrust::make_zip_iterator(thrust::make_tuple(ceventIndex + numEntries, arrayAddress, eventSize)),
            results.begin(),
            *logger);
        host_results = results;

        for(unsigned int j = 0; j < host_results.size(); ++j) {
            values[1 + i].push_back(host_results[j]);
        }
    }

    // restore previous fit control
    setFitControl(fc);

    return values;
}

} // namespace GooFit
