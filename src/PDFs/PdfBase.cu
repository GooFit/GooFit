#include <goofit/BinnedDataSet.h>
#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/Log.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PdfBase.h>
#include <goofit/Version.h>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

// This is code that belongs to the PdfBase class, that is,
// it is common across all implementations. But it calls on device-side
// functions, and due to the nvcc translation-unit limitations, it cannot
// sit in its own object file; it must go in the CUDAglob.cu. So it's
// off on its own in this inline-cuda file, which GooPdf.cu
// should include.

__host__ void PdfBase::copyParams() {
    // Copies values of Variable objects
    std::vector<Variable> pars = getParameters();
    std::vector<double> values;

    updateParameters();
}

__host__ void PdfBase::copyNormFactors() const {
    cudaDeviceSynchronize(); // Ensure normalization integrals are finished
}

__host__ void PdfBase::initializeIndices() {
    // Structure of the individual index array: Number of parameters, then the indices
    // requested by the subclass (which will be interpreted by the subclass kernel),
    // then the number of observables, then the observable indices. Notice that the
    // observable indices are not set until 'setIndices' is called, usually from setData;
    // here we only reserve space for them by setting totalParams.
    // This is to allow index sharing between PDFs - all the PDFs must be constructed
    // before we know what observables exist.

    GOOFIT_DEBUG("PdfBase::initializeIndices for \"{}\"", getName());

    // We are formulating the four buffers.  We first populate the lists based on what is currently in the appropriate
    // PDF buffers.  NOTE: This is not the actual buffer layout!!!!!!!  We are allocating placeholders for each PDF, but
    // they will be re-arranged when the tree is flattened.  Please follow  recursiveSetIndices()

    // stick placeholders into our parameter array
    host_parameters[totalParameters] = parametersList.size();
    parametersIdx                    = totalParameters;
    totalParameters++;
    for(auto &i : parametersList) {
        host_parameters[totalParameters] = i.getValue();
        totalParameters++;
    }

    // stick placeholders into our constants array
    host_constants[totalConstants] = constantsList.size();
    constantsIdx                   = totalConstants;
    totalConstants++;
    for(double i : constantsList) {
        host_constants[totalConstants] = i;
        totalConstants++;
    }

    // stick placeholders into our observable array
    host_observables[totalObservables] = observablesList.size();
    observablesIdx                     = totalObservables;
    totalObservables++;
    for(auto &i : observablesList) {
        host_observables[totalObservables] = i.getValue();
        totalObservables++;
    }

    // stick placeholders into our normalization array
    host_normalizations[totalNormalizations] = 1;
    normalIdx                                = totalNormalizations;
    totalNormalizations++;
    host_normalizations[totalNormalizations] = cachedNormalization;
    totalNormalizations++;

    if(totalParameters >= GOOFIT_MAXPAR)
        throw GooFit::GeneralError("{}: Set too many parameters, GooFit array more than {}. Increase max at compile "
                                   "time with -DGOOFIT_MAXPAR=N.",
                                   getName(),
                                   GOOFIT_MAXPAR);

    // we rely on GooPdf::set to copy these values, this copyies every PDF which is unnecessary.
    // MEMCPY_TO_SYMBOL(paramIndices, host_indices, totalParams*sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::recursiveSetIndices() {
    if(reflex_name_.empty() || function_ptr_ == nullptr)
        throw GeneralError("A PDF must either provide a function name and"
                           " function pointer or override recursiveSetIndices\n"
                           "Called by {}\n"
                           "Make sure initilize is not called before registerFunction!",
                           getName());

    host_fcn_ptr = function_ptr_;

    GOOFIT_DEBUG("host_function_table[{}] = {} for \"{}\" from PDFBase::recursiveSetIndices",
                 num_device_functions,
                 reflex_name_,
                 getName());

    if(num_device_functions >= GOOFIT_MAXFUNC)
        throw GeneralError("Too many device functions! Set GOOFIT_MAXFUNC to a larger value than {}", GOOFIT_MAXFUNC);

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;
    populateArrays();
}

__host__ void PdfBase::updateVariable(Variable var, fptype newValue) {
    for(auto &i : parametersList) {
        if(i.getName() == var.getName())
            i.setValue(newValue);
    }

    for(auto &component : components)
        component->updateVariable(var, newValue);
}

__host__ void PdfBase::updateParameters() {
    for(int i = 0; i < parametersList.size(); i++)
        host_parameters[parametersIdx + i + 1] = parametersList[i].getValue();

    for(auto &component : components)
        component->updateParameters();

    // we need to memcpy to device.
    MEMCPY_TO_SYMBOL(d_parameters, host_parameters, totalParameters * sizeof(fptype), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::populateArrays() {
    // populate all the arrays
    GOOFIT_TRACE("Populating Arrays for {}", getName());

    // reconfigure the host_parameters array with the new indexing scheme.
    GOOFIT_TRACE("host_parameters[{}] = {}", totalParameters, parametersList.size());
    host_parameters[totalParameters] = parametersList.size();
    parametersIdx                    = totalParameters;
    totalParameters++;
    for(auto &i : parametersList) {
        GOOFIT_TRACE("host_parameters[{}] = {}", totalParameters, i.getValue());
        host_parameters[totalParameters] = i.getValue();
        totalParameters++;
    }

    GOOFIT_TRACE("host_constants[{}] = {}", totalConstants, constantsList.size());
    host_constants[totalConstants] = constantsList.size();
    constantsIdx                   = totalConstants;
    totalConstants++;
    for(double i : constantsList) {
        GOOFIT_TRACE("host_constants[{}] = {}", totalConstants, i);
        host_constants[totalConstants] = i;
        totalConstants++;
    }

    GOOFIT_TRACE("host_observables[{}] = {}", totalObservables, observablesList.size());
    host_observables[totalObservables] = observablesList.size();
    observablesIdx                     = totalObservables;
    totalObservables++;
    for(auto &i : observablesList) {
        GOOFIT_TRACE("host_observables[{}] = {}", totalObservables, i.getIndex());
        host_observables[totalObservables] = i.getIndex();
        totalObservables++;
    }

    GOOFIT_TRACE("host_normalizations[{}] = {}", totalNormalizations, 1);
    host_normalizations[totalNormalizations] = 1;
    normalIdx                                = totalNormalizations++;
    GOOFIT_TRACE("host_normalizations[{}] = {}", totalNormalizations, 0);
    host_normalizations[totalNormalizations] = cachedNormalization;
    totalNormalizations++;

    for(auto &component : components)
        component->recursiveSetIndices();

    generateNormRange();
}

__host__ void PdfBase::setIndices() {
    // Flatten the tree by re-running through the whole PDF with everything zero'd.
    totalParameters      = 0;
    totalConstants       = 0;
    totalObservables     = 0;
    totalNormalizations  = 0;
    num_device_functions = 0;

    // set all associated functions parameters, constants, etc.
    recursiveSetIndices();
}

__host__ void PdfBase::setData(DataSet *data) {
    if(dev_event_array) {
        gooFree(dev_event_array);
        cudaDeviceSynchronize();
        dev_event_array  = nullptr;
        m_iEventsPerTask = 0;
    }

    // fine to set static integer value for variables
    // (prod->gauss [somVar->getObservableIndex()]
    //    ->exp [someVar->getObservableIndex()]
    setupObservables();
    setIndices();

    // Tracking the data structure
    data_ = data;

    UnbinnedDataSet *unbinned_data;
    BinnedDataSet *binned_data;

    // Do nothing if passed a nullptr (makes setData(getData()) safe)
    if(data == nullptr) {
        return;
    } else if((unbinned_data = dynamic_cast<UnbinnedDataSet *>(data))) {
        numEntries = data->getNumEvents();
        numEvents  = numEntries;

        size_t dimensions = observablesList.size();

#ifdef GOOFIT_MPI
        // This fetches our rank and the total number of processes in the MPI call
        int myId, numProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);

        int perTask = numEvents / numProcs;

        // This will track for a given rank where they will start and how far they will go
        int *counts        = new int[numProcs];
        int *displacements = new int[numProcs];

        for(int i = 0; i < numProcs - 1; i++)
            counts[i] = perTask;

        counts[numProcs - 1] = numEntries - perTask * (numProcs - 1);

        displacements[0] = 0;

        for(int i = 1; i < numProcs; i++)
            displacements[i] = displacements[i - 1] + counts[i - 1];

#endif

        auto *host_array = new fptype[numEntries * dimensions];

#ifdef GOOFIT_MPI
        // This is an array to track if we need to re-index the observable
        int fixme[observablesList.size()];
        memset(fixme, 0, sizeof(int) * observablesList.size());

        for(int i = 0; i < observables.size(); i++) {
            // We are casting the observable to a CountVariable
            EventNumber *c = dynamic_cast<EventNumber *>(observablesList[i]);

            // if it is true re-index
            if(c)
                fixme[i] = 1;
        }

#endif

        // Transfer into our whole buffer
        for(int i = 0; i < numEntries; ++i) {
            for(const Observable &v : observablesList) {
                fptype currVal                            = unbinned_data->getValue(v, i);
                host_array[i * dimensions + v.getIndex()] = currVal;
            }
        }

#ifdef GOOFIT_MPI

        // We will go through all of the events and re-index if appropriate
        for(int i = 1; i < numProcs; i++) {
            for(int j = 0; j < counts[i]; j++) {
                for(int k = 0; k < dimensions; k++) {
                    if(fixme[k] > 0)
                        host_array[(j + displacements[i]) * dimensions + k] = float(j);
                }
            }
        }

        int mystart = displacements[myId];
        int myend   = mystart + counts[myId];
        int mycount = myend - mystart;

        gooMalloc((void **)&dev_event_array, dimensions * mycount * sizeof(fptype));
        MEMCPY(dev_event_array,
               host_array + mystart * dimensions,
               dimensions * mycount * sizeof(fptype),
               cudaMemcpyHostToDevice);
        // MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
        delete[] host_array;

        setNumPerTask(this, mycount);

        delete[] counts;
        delete[] displacements;
#else
        gooMalloc((void **)&dev_event_array, dimensions * numEntries * sizeof(fptype));
        MEMCPY(dev_event_array, host_array, dimensions * numEntries * sizeof(fptype), cudaMemcpyHostToDevice);
        delete[] host_array;
#endif
        MEMCPY_TO_SYMBOL(c_totalEvents, &numEntries, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
    } else if((binned_data = dynamic_cast<BinnedDataSet *>(data))) {
        numEvents      = 0;
        numEntries     = binned_data->getNumBins();
        int dimensions = 2 + observablesList.size(); // Bin center (x,y, ...), bin value, and bin volume.

        if(!fitControl->binnedFit())
            setFitControl(std::make_shared<BinnedNllFit>());

#ifdef GOOFIT_MPI
        // This fetches our rank and the total number of processes in the MPI call
        int myId, numProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);

        int perTask = numEvents / numProcs;

        // This will track for a given rank where they will start and how far they will go
        int *counts        = new int[numProcs];
        int *displacements = new int[numProcs];

        for(int i = 0; i < numProcs - 1; i++)
            counts[i] = perTask;

        counts[numProcs - 1] = numEntries - perTask * (numProcs - 1);

        displacements[0] = 0;

        for(int i = 1; i < numProcs; i++)
            displacements[i] = displacements[i - 1] + counts[i - 1];
#endif

        auto *host_array = new fptype[numEntries * dimensions];

#ifdef GOOFIT_MPI
        // This is an array to track if we need to re-index the observable
        int fixme[observables.size()];
        memset(fixme, 0, sizeof(int) * observablesList.size());

        for(int i = 0; i < observablesList.size(); i++) {
            // We are casting the observable to a CountVariable
            EventNumber *c = dynamic_cast<EventNumber *>(observablesList[i]);

            // if it is true re-index
            if(c)
                fixme[i] = 1;
        }
#endif

        for(unsigned int i = 0; i < numEntries; ++i) {
            for(const Observable &v : observablesList) {
                host_array[i * dimensions + v.getIndex()] = binned_data->getBinCenter(v, i);
            }

            host_array[i * dimensions + observablesList.size() + 0] = binned_data->getBinContent(i);
            host_array[i * dimensions + observablesList.size() + 1]
                = fitControl->binErrors() ? binned_data->getBinError(i) : binned_data->getBinVolume(i);
            numEvents += binned_data->getBinContent(i);
        }

#ifdef GOOFIT_MPI
        // We will go through all of the events and re-index if appropriate
        for(int i = 1; i < numProcs; i++) {
            for(int j = 0; j < counts[j]; j++) {
                for(int k = 0; k < dimensions; k++) {
                    if(fixme[k] > 0)
                        host_array[(j + displacements[i]) * dimensions + k] = float(j);
                }
            }
        }

        int mystart = displacements[myId];
        int myend   = mystart + counts[myId];
        int mycount = myend - mystart;

        gooMalloc((void **)&dev_event_array, dimensions * mycount * sizeof(fptype));
        MEMCPY(dev_event_array,
               host_array + mystart * dimensions,
               dimensions * mycount * sizeof(fptype),
               cudaMemcpyHostToDevice);
        delete[] host_array;

        setNumPerTask(this, mycount);

        delete[] counts;
        delete[] displacements;
#else
        gooMalloc((void **)&dev_event_array, dimensions * numEntries * sizeof(fptype));
        MEMCPY(dev_event_array, host_array, dimensions * numEntries * sizeof(fptype), cudaMemcpyHostToDevice);
        delete[] host_array;
#endif
        MEMCPY_TO_SYMBOL(c_totalEvents, &numEntries, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
    } else
        throw GooFit::GeneralError("Dataset must be binned or unbinned!");
}

__host__ void PdfBase::generateNormRange() {
    if(normRanges)
        gooFree(normRanges);

    gooMalloc(reinterpret_cast<void **>(&normRanges), 3 * observablesList.size() * sizeof(fptype));

    auto *host_norms = new fptype[3 * observablesList.size()];
    int counter      = 0; // Don't use index in this case to allow for, eg,

    // a single observable whose index is 1; or two observables with indices
    // 0 and 2. Make one array per functor, as opposed to variable, to make
    // it easy to pass MetricTaker a range without worrying about which parts
    // to use.
    for(Observable &v : observablesList) {
        host_norms[3 * counter + 0] = v.getLowerLimit();
        host_norms[3 * counter + 1] = v.getUpperLimit();
        host_norms[3 * counter + 2] = integrationBins > 0 ? integrationBins : v.getNumBins();
        GOOFIT_TRACE("host_norms[{}] = {}", 3 * counter + 0, host_norms[3 * counter + 0]);
        GOOFIT_TRACE("host_norms[{}] = {}", 3 * counter + 1, host_norms[3 * counter + 1]);
        GOOFIT_TRACE("host_norms[{}] = {}", 3 * counter + 2, host_norms[3 * counter + 2]);
        counter++;
    }

    MEMCPY(normRanges, host_norms, 3 * observablesList.size() * sizeof(fptype), cudaMemcpyHostToDevice);
    delete[] host_norms;
}

void PdfBase::clearCurrentFit() {
    totalParameters = 0;
    gooFree(dev_event_array);
    dev_event_array = nullptr;
}

cudaError_t gooMalloc(void **target, size_t bytes) {
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
    target[0] = malloc(bytes);

    if(target[0])
        return cudaSuccess;
    else
        return cudaErrorMemoryAllocation;

#else
    return cudaMalloc(target, bytes);
#endif
}

cudaError_t gooFree(void *ptr) {
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
    free(ptr);
    return cudaSuccess;
#else
    return cudaFree(ptr);
#endif
}
} // namespace GooFit
