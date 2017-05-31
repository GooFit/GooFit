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

#ifdef ROOT_FOUND
#include <TH1D.h>
#endif

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

// These variables are either function-pointer related (thus specific to this implementation)
// or constrained to be in the CUDAglob translation unit by nvcc limitations; otherwise they
// would be in PdfBase.

// Device-side, translation-unit constrained.

__constant__ fptype cudaArray[maxParams];
__constant__ unsigned int paramIndices[maxParams];
__constant__ fptype functorConstants[maxParams];
__constant__ fptype normalisationFactors[maxParams];

// For debugging

__constant__ int callnumber;
__constant__ int gpuDebug;
__constant__ unsigned int debugParamIndex;
__device__ int internalDebug1 = -1;
__device__ int internalDebug2 = -1;
__device__ int internalDebug3 = -1;
int cpuDebug                  = 0;

#ifdef PROFILING
__device__ fptype timeHistogram[10000];
fptype host_timeHist[10000];
#endif

// Function-pointer related.
__device__ void *device_function_table[200];
// Not clear why this cannot be __constant__, but it causes crashes to declare it so.

void *host_function_table[200];
unsigned int num_device_functions = 0;
std::map<void *, int> functionAddressToDeviceIndexMap;

// For use in debugging memory issues
void printMemoryStatus(std::string file, int line) {
    size_t memfree  = 0;
    size_t memtotal = 0;
    cudaDeviceSynchronize();

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    cudaMemGetInfo(&memfree, &memtotal);
#endif
    cudaDeviceSynchronize();
    std::cout << "Memory status " << file << " " << line << " Free " << memfree << " Total " << memtotal << " Used "
              << (memtotal - memfree) << std::endl;
}

__device__ fptype calculateEval(fptype rawPdf, fptype *evtVal, unsigned int par) {
    // Just return the raw PDF value, for use in (eg) normalisation.
    return rawPdf;
}

__device__ fptype calculateNLL(fptype rawPdf, fptype *evtVal, unsigned int par) {
    // if ((10 > callnumber) && (THREADIDX < 10) && (BLOCKIDX == 0)) cuPrintf("calculateNll %i %f %f %f\n", callnumber,
    // rawPdf, normalisationFactors[par], rawPdf*normalisationFactors[par]);
    // if (THREADIDX < 50) printf("Thread %i %f %f\n", THREADIDX, rawPdf, normalisationFactors[par]);
    rawPdf *= normalisationFactors[par];
    return rawPdf > 0 ? -log(rawPdf) : 0;
}

__device__ fptype calculateProb(fptype rawPdf, fptype *evtVal, unsigned int par) {
    // Return probability, ie normalized PDF value.
    return rawPdf * normalisationFactors[par];
}

__device__ fptype calculateBinAvg(fptype rawPdf, fptype *evtVal, unsigned int par) {
    rawPdf *= normalisationFactors[par];
    rawPdf *= evtVal[1]; // Bin volume

    // Log-likelihood of numEvents with expectation of exp is (-exp + numEvents*ln(exp) - ln(numEvents!)).
    // The last is constant, so we drop it; and then multiply by minus one to get the negative log-likelihood.
    if(rawPdf > 0) {
        fptype expEvents = functorConstants[0] * rawPdf;
        return (expEvents - evtVal[0] * log(expEvents));
    }

    return 0;
}

__device__ fptype calculateBinWithError(fptype rawPdf, fptype *evtVal, unsigned int par) {
    // In this case interpret the rawPdf as just a number, not a number of events.
    // Do not divide by integral over phase space, do not multiply by bin volume,
    // and do not collect 200 dollars. evtVal should have the structure (bin entry, bin error).
    // printf("[%i, %i] ((%f - %f) / %f)^2 = %f\n", BLOCKIDX, THREADIDX, rawPdf, evtVal[0], evtVal[1], pow((rawPdf -
    // evtVal[0]) / evtVal[1], 2));
    rawPdf -= evtVal[0]; // Subtract observed value.
    rawPdf /= evtVal[1]; // Divide by error.
    rawPdf *= rawPdf;
    return rawPdf;
}

__device__ fptype calculateChisq(fptype rawPdf, fptype *evtVal, unsigned int par) {
    rawPdf *= normalisationFactors[par];
    rawPdf *= evtVal[1]; // Bin volume

    return POW2(rawPdf * functorConstants[0] - evtVal[0]) / (evtVal[0] > 1 ? evtVal[0] : 1);
}

__device__ device_metric_ptr ptr_to_Eval         = calculateEval;
__device__ device_metric_ptr ptr_to_NLL          = calculateNLL;
__device__ device_metric_ptr ptr_to_Prob         = calculateProb;
__device__ device_metric_ptr ptr_to_BinAvg       = calculateBinAvg;
__device__ device_metric_ptr ptr_to_BinWithError = calculateBinWithError;
__device__ device_metric_ptr ptr_to_Chisq        = calculateChisq;

void *host_fcn_ptr = nullptr;

void *getMetricPointer(std::string name) {
#define CHOOSE_PTR(ptrname)                                                                                            \
    if(name == #ptrname)                                                                                               \
        GET_FUNCTION_ADDR(ptrname);
    host_fcn_ptr = nullptr;
    CHOOSE_PTR(ptr_to_Eval);
    CHOOSE_PTR(ptr_to_NLL);
    CHOOSE_PTR(ptr_to_Prob);
    CHOOSE_PTR(ptr_to_BinAvg);
    CHOOSE_PTR(ptr_to_BinWithError);
    CHOOSE_PTR(ptr_to_Chisq);

    if(host_fcn_ptr == nullptr)
        throw GooFit::GeneralError("host_fcn_ptr is nullptr");

    return host_fcn_ptr;
#undef CHOOSE_PTR
}

GooPdf::GooPdf(Variable *x, std::string n)
    : PdfBase(x, n)
    , logger(nullptr) {
    // std::cout << "Created " << n << std::endl;
}

__host__ int GooPdf::findFunctionIdx(void *dev_functionPtr) {
    // Code specific to function-pointer implementation
    auto localPos = functionAddressToDeviceIndexMap.find(dev_functionPtr);

    if(localPos != functionAddressToDeviceIndexMap.end()) {
        return (*localPos).second;
    }

    int fIdx                                         = num_device_functions;
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

__host__ void GooPdf::initialize(std::vector<unsigned int> pindices, void *dev_functionPtr) {
    if(!fitControl)
        setFitControl(new UnbinnedNllFit());

    // MetricTaker must be created after PdfBase initialisation is done.
    PdfBase::initializeIndices(pindices);

    functionIdx = findFunctionIdx(dev_functionPtr);
    setMetrics();
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

__host__ void GooPdf::setMetrics() {
    if(logger)
        delete logger;

    logger = new MetricTaker(this, getMetricPointer(fitControl->getMetric()));
}

__host__ double GooPdf::sumOfNll(int numVars) const {
    static thrust::plus<double> cudaPlus;
    thrust::constant_iterator<int> eventSize(numVars);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    double dummy = 0;

    // if (host_callnumber >= 2) GooFit::abort(__FILE__, __LINE__, getName() + " debug abort", this);
    thrust::counting_iterator<int> eventIndex(0);

    double ret;
#ifdef GOOFIT_MPI
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    goofit_policy my_policy;
    double r = thrust::transform_reduce(
        my_policy,
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#else
    double r = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#endif

    MPI_Allreduce(&r, &ret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    goofit_policy my_policy;
    ret = thrust::transform_reduce(
        my_policy,
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#else
    ret = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
        *logger,
        dummy,
        cudaPlus);
#endif
#endif
    return ret;
}

__host__ double GooPdf::calculateNLL() const {
    // if (cpuDebug & 1) std::cout << getName() << " entering calculateNLL (" << host_callnumber << ")" << std::endl;

    // MEMCPY_TO_SYMBOL(callnumber, &host_callnumber, sizeof(int));
    // int oldMask = cpuDebug;
    // if (0 == host_callnumber) setDebugMask(0, false);
    // std::cout << "Start norm " << getName() << std::endl;
    normalize();
    // std::cout << "Norm done\n";
    // if ((0 == host_callnumber) && (1 == oldMask)) setDebugMask(1, false);

    // if (cpuDebug & 1) {
    // std::cout << "Norm factors: ";
    // for (int i = 0; i < totalParams; ++i) std::cout << host_normalisation[i] << " ";
    // std::cout << std::endl;
    //}

    if(host_normalisation[parameters] <= 0)
        GooFit::abort(__FILE__, __LINE__, getName() + " non-positive normalisation", this);

    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize(); // Ensure normalisation integrals are finished

    int numVars = observables.size();

    if(fitControl->binnedFit()) {
        numVars += 2;
        numVars *= -1;
    }

    fptype ret = sumOfNll(numVars);

    if(0 == ret)
        GooFit::abort(__FILE__, __LINE__, getName() + " zero NLL", this);

    // if (cpuDebug & 1) std::cout << "Full NLL " << host_callnumber << " : " << 2*ret << std::endl;
    // setDebugMask(0);

    // if ((cpuDebug & 1) && (host_callnumber >= 1)) GooFit::abort(__FILE__, __LINE__, getName() + " debug abort",
    // this);
    return 2 * ret;
}

__host__ std::vector<fptype> GooPdf::evaluateAtPoints(Variable *var) {
    copyParams();
    normalize();
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);
    UnbinnedDataSet tempdata(observables);

    double step = var->getBinSize();

    for(int i = 0; i < var->getNumBins(); ++i) {
        var->setValue(var->getLowerLimit() + (i + 0.5) * step);
        tempdata.addEvent();
    }

    setData(&tempdata);

    thrust::counting_iterator<int> eventIndex(0);
    thrust::constant_iterator<int> eventSize(observables.size());
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::device_vector<fptype> results(var->getNumBins());

    MetricTaker evalor(this, getMetricPointer("ptr_to_Eval"));
#ifdef GOOFIT_MPI
    thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
        results.begin(),
        evalor);
#else
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                      results.begin(),
                      evalor);
#endif

    // Note, This is not fully realized with MPI.  We need to copy each 'results' buffer to each other 'MPI_Scatterv',
    // then we can do the rest.
    thrust::host_vector<fptype> h_results = results;
    std::vector<fptype> res;
    res.resize(var->getNumBins());

    for(int i = 0; i < var->getNumBins(); ++i) {
        res[i] = h_results[i] * host_normalisation[parameters];
    }

    return res;
}

__host__ void GooPdf::scan(Variable *var, std::vector<fptype> &values) {
    fptype step = var->getUpperLimit();
    step -= var->getLowerLimit();
    step /= var->getNumBins();
    values.clear();

    for(fptype v = var->getLowerLimit() + 0.5 * step; v < var->getUpperLimit(); v += step) {
        var->setValue(v);
        copyParams();
        fptype curr = calculateNLL();
        values.push_back(curr);
    }
}

// TODO: is this needed?
__host__ void GooPdf::setParameterConstantness(bool constant) {
    std::vector<Variable *> pars = getParameters();

    for(Variable *p : pars) {
        p->setFixed(constant);
    }
}

__host__ fptype GooPdf::getValue() {
    // Returns the value of the PDF at a single point.
    // Execute redundantly in all threads for OpenMP multiGPU case
    copyParams();
    normalize();
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    UnbinnedDataSet point(observables);
    point.addEvent();
    setData(&point);

    thrust::counting_iterator<int> eventIndex(0);
    thrust::constant_iterator<int> eventSize(observables.size());
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::device_vector<fptype> results(1);

    MetricTaker evalor(this, getMetricPointer("ptr_to_Eval"));
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + 1, arrayAddress, eventSize)),
                      results.begin(),
                      evalor);

    return results[0];
}

__host__ fptype GooPdf::normalize() const {
    if(!fitControl->metricIsPdf()) {
        GOOFIT_TRACE("{}: metricIsPdf, returning 1", getName());
        host_normalisation[parameters] = 1.0;
        return 1.0;
    }

    fptype ret = 1;

    if(hasAnalyticIntegral()) {
        // Loop goes only over observables of this PDF.
        for(Variable *v : observables) {
            GOOFIT_TRACE("{}: Analytically integrating over {}", getName(), v->getName());
            ret *= integrate(v->getLowerLimit(), v->getUpperLimit());
        }

        host_normalisation[parameters] = 1.0 / ret;
        GOOFIT_TRACE("{}: Param {} integral is = {}", getName(), parameters, ret);

        return ret;
    }

    GOOFIT_TRACE("{}, Computing integral without analytic help", getName());

    int totalBins = 1;

    for(Variable *v : observables) {
        ret *= v->getUpperLimit() - v->getLowerLimit();
        totalBins *= integrationBins > 0 ? integrationBins : v->getNumBins();

        GOOFIT_TRACE("Total bins {} due to {} {} {}", totalBins, v->getName(), integrationBins, v->getNumBins());
    }

    ret /= totalBins;

    fptype dummy = 0;
    static thrust::plus<fptype> cudaPlus;
    thrust::constant_iterator<fptype *> arrayAddress(normRanges);
    thrust::constant_iterator<int> eventSize(observables.size());
    thrust::counting_iterator<int> binIndex(0);

    fptype sum;
#ifdef GOOFIT_MPI
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    goofit_policy my_policy;
    fptype s = thrust::transform_reduce(
        my_policy,
        thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
        thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
        *logger,
        dummy,
        cudaPlus);
#else
    fptype s = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
        thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
        *logger,
        dummy,
        cudaPlus);
#endif

    MPI_Allreduce(&s, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    goofit_policy my_policy;
    sum = thrust::transform_reduce(
        my_policy,
        thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
        thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
        *logger,
        dummy,
        cudaPlus);
#else
    sum = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
        thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
        *logger,
        dummy,
        cudaPlus);
#endif
#endif

    if(std::isnan(sum)) {
        GooFit::abort(__FILE__, __LINE__, getName() + " NaN in normalisation", this);
    } else if(0 >= sum) {
        GooFit::abort(__FILE__, __LINE__, "Non-positive normalisation", this);
    }

    ret *= sum;

    if(0 == ret)
        GooFit::abort(__FILE__, __LINE__, "Zero integral");

    GOOFIT_TRACE("{}: Param {} integral is ~= {}", getName(), parameters, ret);
    host_normalisation[parameters] = 1.0 / ret;
    return ret;
}

#ifdef PROFILING
__constant__ fptype conversion = (1.0 / CLOCKS_PER_SEC);
__device__ fptype callFunction(fptype *eventAddress, unsigned int functionIdx, unsigned int paramIdx) {
    clock_t start = clock();
    fptype ret    = (*(reinterpret_cast<device_function_ptr>(device_function_table[functionIdx])))(
        eventAddress, cudaArray, paramIndices + paramIdx);
    clock_t stop = clock();

    if((0 == THREADIDX + BLOCKIDX) && (stop > start)) {
        // Avoid issue when stop overflows and start doesn't.
        timeHistogram[functionIdx * 100 + paramIdx] += ((stop - start) * conversion);
        // printf("Clock: %li %li %li | %u %f\n", (long) start, (long) stop, (long) (stop - start), functionIdx,
        // timeHistogram[functionIdx]);
    }

    return ret;
}
#else
__device__ fptype callFunction(fptype *eventAddress, unsigned int functionIdx, unsigned int paramIdx) {
    return (*(reinterpret_cast<device_function_ptr>(device_function_table[functionIdx])))(
        eventAddress, cudaArray, paramIndices + paramIdx);
}
#endif

__host__ std::vector<std::vector<fptype>> GooPdf::getCompProbsAtDataPoints() {
    copyParams();
    // double overall =
    normalize();
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    int numVars = observables.size();

    if(fitControl->binnedFit()) {
        numVars += 2;
        numVars *= -1;
    }

    thrust::device_vector<fptype> results(numEntries);
    thrust::constant_iterator<int> eventSize(numVars);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);
    MetricTaker evalor(this, getMetricPointer("ptr_to_Prob"));
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                      results.begin(),
                      evalor);
    std::vector<std::vector<fptype>> values;
    values.resize(components.size() + 1);
    thrust::host_vector<fptype> host_results = results;

    for(unsigned int i = 0; i < host_results.size(); ++i) {
        values[0].push_back(host_results[i]);
    }

    for(unsigned int i = 0; i < components.size(); ++i) {
        MetricTaker compevalor(components[i], getMetricPointer("ptr_to_Prob"));
        thrust::counting_iterator<int> ceventIndex(0);
        thrust::transform(
            thrust::make_zip_iterator(thrust::make_tuple(ceventIndex, arrayAddress, eventSize)),
            thrust::make_zip_iterator(thrust::make_tuple(ceventIndex + numEntries, arrayAddress, eventSize)),
            results.begin(),
            compevalor);
        host_results = results;

        for(unsigned int j = 0; j < host_results.size(); ++j) {
            values[1 + i].push_back(host_results[j]);
        }
    }
    return values;
}

// Utility function to make a grid of any dimisinion
__host__ void make_a_grid(std::vector<Variable *> ret, UnbinnedDataSet &grid) {
    if(ret.empty()) {
        grid.addEvent();
        return;
    }

    Variable *var = ret.back();
    ret.pop_back(); // safe because this is a copy

    for(int i = 0; i < var->getNumBins(); ++i) {
        double step = (var->getUpperLimit() - var->getLowerLimit()) / var->getNumBins();
        var->setValue(var->getLowerLimit() + (i + 0.5) * step);
        make_a_grid(ret, grid);
    }
}

__host__ UnbinnedDataSet GooPdf::makeGrid() {
    std::vector<Variable *> ret = getObservables();

    UnbinnedDataSet grid{ret};

    make_a_grid(ret, grid);

    return grid;
}

// still need to add OpenMP/multi-GPU code here
__host__ void GooPdf::transformGrid(fptype *host_output) {
    generateNormRange();
    // normalize();
    int totalBins = 1;

    for(Variable *v : observables) {
        totalBins *= v->getNumBins();
    }

    thrust::constant_iterator<fptype *> arrayAddress(normRanges);
    thrust::constant_iterator<int> eventSize(observables.size());
    thrust::counting_iterator<int> binIndex(0);
    thrust::device_vector<fptype> d_vec;
    d_vec.resize(totalBins);

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
                      thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, eventSize, arrayAddress)),
                      d_vec.begin(),
                      *logger);

    thrust::host_vector<fptype> h_vec = d_vec;

    for(unsigned int i = 0; i < totalBins; ++i)
        host_output[i] = h_vec[i];
}

__host__ void GooPdf::setFitControl(FitControl *const fc, bool takeOwnerShip) {
    for(auto &component : components) {
        component->setFitControl(fc, false);
    }

    if((fitControl) && (fitControl->getOwner() == this)) {
        delete fitControl;
    }

    fitControl = fc;

    if(takeOwnerShip) {
        fitControl->setOwner(this);
    }

    setMetrics();
}

#ifdef ROOT_FOUND
__host__ TH1D *GooPdf::plotToROOT(Variable *var, double normFactor, std::string name) {
    if(name.empty())
        name = getName() + "_hist";

    auto ret = new TH1D(name.c_str(), "", var->getNumBins(), var->getLowerLimit(), var->getUpperLimit());
    std::vector<fptype> binValues = evaluateAtPoints(var);

    double pdf_int = 0;

    for(int i = 0; i < var->getNumBins(); ++i) {
        pdf_int += binValues[i];
    }

    for(int i = 0; i < var->getNumBins(); ++i)
        ret->SetBinContent(i + 1, binValues[i] * normFactor / pdf_int / var->getBinSize());
    return ret;
}
#endif
} // namespace GooFit
