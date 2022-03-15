#include <goofit/BinnedDataSet.h>
#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/Log.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>
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

__host__ void PdfBase::pre_run() {
    GOOFIT_TRACE("GooPdf::pre_run");
    host_function_table.sync(d_function_table);
    host_parameters.sync(d_parameters);
    host_constants.sync(d_constants);
    host_observables.sync(d_observables);
    host_normalizations.sync(d_normalizations);
}

__host__ void PdfBase::pre_call() {
    GOOFIT_TRACE("GooPdf::pre_call");
    host_parameters.smart_sync(d_parameters);
}

__host__ void PdfBase::copyParams() {
    // Copies values of Variable objects
    std::vector<Variable> pars = getParameters();
    std::vector<double> values;

    updateParameters();
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
    parametersIdx = host_parameters.size();
    host_parameters.push_back(parametersList.size());
    host_parameter_name.push_back("nParam");
    for(auto &i : parametersList) {
        host_parameters.push_back(i.getValue());
        host_parameter_name.push_back(i.getName());
    }

    // stick placeholders into our constants array
    constantsIdx = host_constants.size();
    host_constants.push_back(constantsList.size());
    for(double i : constantsList) {
        host_constants.push_back(i);
    }

    // stick placeholders into our observable array
    observablesIdx = host_observables.size();
    host_observables.push_back(observablesList.size());
    host_observable_name.push_back("nObs");
    for(auto &i : observablesList) {
        host_observables.push_back(i.getValue());
        host_observable_name.push_back(i.getName());
    }

    // stick placeholders into our normalization array
    normalIdx = host_normalizations.size();
    host_normalizations.push_back(1.0);
    host_normalizations.push_back(0.0);

    // we rely on GooPdf::set to copy these values, this copyies every PDF which is unnecessary.
    // MEMCPY_TO_SYMBOL(paramIndices, host_indices, totalParams*sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::recursiveSetIndices() {
    if(reflex_name_.empty() || function_ptr_ == nullptr)
        throw GeneralError("A PDF must either provide a function name and"
                           " function pointer or override recursiveSetIndices\n"
                           "Called by {}\n"
                           "Make sure initialize is not called before registerFunction!",
                           getName());

    host_fcn_ptr = function_ptr_;

    GOOFIT_DEBUG("host_function_table[{}] = {} for \"{}\" from PDFBase::recursiveSetIndices",
                 host_function_table.size(),
                 reflex_name_,
                 getName());

    functionIdx = host_function_table.size();

    host_function_table.push_back(host_fcn_ptr);
    std::string pdfName = getName();
    host_function_name.push_back(pdfName);

    // TODO: This one doesn't have a sync call. Fix and simplify all sync calls!
    populateArrays();
    // mds    PdfBase::status(" in recursiveSetIndices, after populateArrays() ");
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
        host_parameters.at(parametersIdx + i + 1) = parametersList[i].getValue();

    for(auto &component : components)
        component->updateParameters();

    // we need to memcpy to device.
    pre_call();
}

__host__ void PdfBase::status() {
    std::cout << "  ** entered PdfBase::status()  \n";
    for(auto &i : parametersList) {
        auto pName = i.getName();
        std::cout << " parameter name = " << pName << "\n";
    }

    std::cout << "about to print components  \n";
    for(auto &component : components) {
        std::cout << " component = " << component << ",  getName() =  " << getName()
                  << ",  reflex_name_ = " << reflex_name_ << "\n";
    }
    std::cout << "   done printing components  \n";

    auto n_host_function_table = host_function_table.size();
    auto n_host_parameters     = host_parameters.size();
    auto n_host_constants      = host_constants.size();
    auto n_host_observables    = host_observables.size();
    std::cout << " --> n_ host_function_table, parameters, constants, observables = \n"
              << " -->   " << n_host_function_table << "  " << n_host_parameters << "  " << n_host_constants << "  "
              << n_host_observables << "\n\n";

    for(int ii = 0; ii < n_host_function_table; ii++) {
        auto host_function = host_function_table[ii];
        std::cout << "  host_function  " << ii << "  =  " << host_function << "\n";
        auto device_fcn_ptr = d_function_table[ii];
        std::cout << "   device_fcn_ptr       " << device_fcn_ptr << "\n";
        //            auto fIdx = GooPdf::lookUpFunctionIdx(device_fcn_ptr);
        //            std::cout << "   fIdx = " << fIdx << "\n";
    }

    //  host_function->getName() should get the "name"
    //  host_function->reflex_name_

    for(int ii = 0; ii < n_host_parameters; ii++) {
        auto host_parameter = host_parameters[ii];
        std::cout << "  host_parameter  " << ii << "  =  " << host_parameter << "\n";
    }

    for(int ii = 0; ii < n_host_constants; ii++) {
        auto host_constant = host_constants[ii];
        std::cout << "  host_constant  " << ii << "  =  " << host_constant << "\n";
    }

    std::cout << "\n"
              << "  functionPtrToNameMap contains "
              << "\n";
    std::map<void *, std::string>::iterator it = functionPtrToNameMap.begin();
    while(it != functionPtrToNameMap.end()) {
        std::cout << it->first << " :: " << it->second << std::endl;
        it++;
    }

    std::cout << "** about to depart PdfBase::status()  \n";
}

// add a version which takes a string as an argument so we can
// more easily track the program flow
__host__ void PdfBase::status(std::string caller) {
    std::cout << "  ** entered PdfBase::status()  " + caller + "  \n";

    std::cout << "parameter names  \n";
    for(auto &i : parametersList) {
        auto pName  = i.getName();
        auto pValue = i.getValue();
        std::cout << "   parameter name = " << pName << "  =  " << pValue << "\n";
    }

    std::cout << "observable names:  \n";
    for(auto &i : observablesList) {
        auto pName  = i.getName();
        auto pValue = i.getValue();
        std::cout << "   observable name = " << pName << "  =  " << pValue << "\n";
    }

    std::cout << "constants  \n";
    for(auto &cValue : constantsList) {
        std::cout << "   constant value = " << cValue << "\n";
    }

    std::cout << "about to print components  \n";
    for(auto &component : components) {
        std::cout << " component = " << component << ",  getName() =  " << getName()
                  << ",  reflex_name_ = " << reflex_name_ << "\n";
        std::cout << "     *component = " << *component << "\n";
    }
    std::cout << "   done printing components  \n";

    std::cout << "  in PdfBase::status(caller) with caller =  " << caller << "   \n";
    auto size_of_host_function_name = host_function_name.size();
    for(int ii = 0; ii < size_of_host_function_name; ii++) {
        auto aString        = host_function_name[ii];
        auto device_fcn_ptr = d_function_table[ii];

        std::string device_fcn_ptr_name            = "not assigned yet";
        std::map<void *, std::string>::iterator it = functionPtrToNameMap.begin();
        while(it != functionPtrToNameMap.end()) {
            if(device_fcn_ptr == it->first) {
                device_fcn_ptr_name = it->second;
            }
            it++;
        }
        //            auto fIdx = GooPdf::lookUpFunctionIdx(device_fcn_ptr);
        std::cout << "  host_function_name  " << ii << "  =  " << aString << ",  ptr =  " << host_function_table[ii]
                  << "     with device_fcn_ptr_name   " << device_fcn_ptr_name << "\n";
    }
    std::cout << "  \n";
    auto size_of_host_parameters = host_parameters.size();
    for(int ii = 0; ii < size_of_host_parameters; ii++) {
        auto aString = host_parameter_name[ii];
        std::cout << "  host_parameter_name  " << ii << "  =  " << aString << "  = " << host_parameters[ii] << "\n";
    }

    std::cout << "\n"
              << "  functionPtrToNameMap contains "
              << "\n";
    std::map<void *, std::string>::iterator it = functionPtrToNameMap.begin();
    while(it != functionPtrToNameMap.end()) {
        std::cout << it->first << " :: " << it->second << std::endl;
        it++;
    }

    std::cout << " -------- about to exit   PdfBase::status()  " + caller + "  \n";
    std::cout << "\n";
}

__host__ void PdfBase::populateArrays() {
    // populate all the arrays
    GOOFIT_TRACE("Populating Arrays for {}", getName());

    // reconfigure the host_parameters array with the new indexing scheme.
    GOOFIT_TRACE("host_parameters[{}] = {}", host_parameters.size(), parametersList.size());

    parametersIdx = host_parameters.size();
    host_parameters.push_back(parametersList.size());
    host_parameter_name.push_back("nParam");
    for(auto &i : parametersList) {
        GOOFIT_TRACE("host_parameters[{}] = {}", host_parameters.size(), i.getValue());
        host_parameters.push_back(i.getValue());
        host_parameter_name.push_back(i.getName());
    }

    constantsIdx = host_constants.size();
    host_constants.push_back(constantsList.size());

    for(auto &i : constantsList) {
        GOOFIT_TRACE("host_constants[{}] = {}", host_constants.size(), i);
        host_constants.push_back(i);
    }

    observablesIdx = host_observables.size();
    host_observables.push_back(observablesList.size());
    host_observable_name.push_back("nObs");
    for(auto &i : observablesList) {
        GOOFIT_TRACE("host_observables[{}] = {}", host_observables.size(), i.getValue());
        host_observables.push_back(i.getIndex());
        host_observable_name.push_back(i.getName());
    }

    GOOFIT_TRACE("host_normalizations[{}] = {}", host_normalizations.size(), 1);
    normalIdx = host_normalizations.size();
    host_normalizations.push_back(1.0);
    GOOFIT_TRACE("host_normalizations[{}] = {}", host_normalizations.size(), 0);
    host_normalizations.push_back(cachedNormalization);

    for(auto &component : components)
        component->recursiveSetIndices();

    generateNormRange();
}

__host__ void PdfBase::setIndices() {
    // Flatten the tree by re-running through the whole PDF with everything zero'd.
    host_parameters.clear();
    host_constants.clear();
    host_observables.clear();
    host_normalizations.clear();
    host_function_table.clear();

    host_function_name.clear();
    host_parameter_name.clear();
    host_observable_name.clear();

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

        int idx = 0;
        for(auto o : observablesList) {
            // We are casting the observable to a EventNumber
            fixme[idx] = o.isEventNumber();
            idx++;
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
        int fixme[observablesList.size()];
        memset(fixme, 0, sizeof(int) * observablesList.size());

        int idx = 0;
        for(auto o : observablesList) {
            // if it is true re-index
            fixme[idx] = o.isEventNumber();
            idx++;
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

auto gooMalloc(void **target, size_t bytes) -> cudaError_t {
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

auto gooFree(void *ptr) -> cudaError_t {
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
    free(ptr);
    return cudaSuccess;
#else
    return cudaFree(ptr);
#endif
}
} // namespace GooFit
