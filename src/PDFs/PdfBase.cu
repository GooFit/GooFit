#include "goofit/GlobalCudaDefines.h"
#include "goofit/PdfBase.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Error.h"
#include "goofit/Log.h"
#include "goofit/BinnedDataSet.h"
#include "goofit/FitControl.h"

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

__host__ void PdfBase::copyParams(const std::vector<double>& pars) const {
    // copyParams method performs eponymous action!

    //for(unsigned int i = 0; i < pars.size(); ++i) {
    //    host_parameters[parameterIdx + i + 1] = pars[i];

    //    if(std::isnan(host_parameters[i])) {
    //        std::cout << " agh, parameter is NaN, die " << i << std::endl;
    //        GooFit::abort(__FILE__, __LINE__, "NaN in parameter");
    //    }
    //}

    for (int i = 0; i < parametersList.size (); i++)
    {
        GOOFIT_TRACE("fitter index {}", parametersList[i]->getFitterIndex());
        host_parameters[parametersIdx + i + 1] = pars[parametersList[i]->getFitterIndex()];
    }

    MEMCPY(d_parameters, host_parameters, totalParameters*sizeof(fptype), cudaMemcpyHostToDevice);

    //recursiveSetIndices ();
    //MEMCPY_TO_SYMBOL(cudaArray, host_params, pars.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::copyParams() {
    // Copies values of Variable objects
    std::vector<Variable*> pars = getParameters();
    std::vector<double> values;

    //for(Variable* v : pars) {
    //    int index = v->getIndex();

    //    if(index >= (int) values.size())
    //        values.resize(index + 1);

    //    values[index] = v->getValue();
    //}

    copyParams(values);

    updateParameters ();
}

__host__ void PdfBase::copyNormFactors() const {
    //MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize(); // Ensure normalisation integrals are finished
}

__host__ void PdfBase::initialiseIndices(std::vector<unsigned int> pindices) {
    // Structure of the individual index array: Number of parameters, then the indices
    // requested by the subclass (which will be interpreted by the subclass kernel),
    // then the number of observables, then the observable indices. Notice that the
    // observable indices are not set until 'setIndices' is called, usually from setData;
    // here we only reserve space for them by setting totalParams.
    // This is to allow index sharing between PDFs - all the PDFs must be constructed
    // before we know what observables exist.
    
    GOOFIT_DEBUG("Adding space for {} indices for {}", pindices.size(), getName());

    //TODO:We need to check parameters, constants, observables, and norms
    if(totalParameters + pindices.size() >= maxParams)
        throw GooFit::GeneralError("totalParams {} + pindices {} must be less than {}", totalParameters, pindices.size(), maxParams);

    //We are formulating the four buffers.  We first populate the lists based on what is currently in the appropriate PDF buffers.
    //NOTE: This is not the actual buffer layout!!!!!!!
    //We are allocating placeholders for each PDF, but they will be re-arranged when the tree is flattened.  Please follow
    //recursiveSetIndices()

    //stick placeholders into our parameter array 
    host_parameters[totalParameters] = parametersList.size();
    parametersIdx = totalParameters; totalParameters++; 
    for (int i = 0 ; i < parametersList.size (); i++) {
        host_parameters[totalParameters] = parametersList[i]->getValue (); totalParameters++;
    }

    //stick placeholders into our observable array
    host_observables[totalObservables] = observablesList.size ();
    observablesIdx = totalObservables; totalObservables++;
    for (int i = 0 ; i < observablesList.size (); i++) {
        host_observables[totalObservables] = observablesList[i]->getValue (); totalObservables++;
    }

    //stick placeholders into our constants array
    host_constants[totalConstants] = constantsList.size ();
    constantsIdx = totalConstants; totalConstants++;
    for (int i = 0 ; i < constantsList.size (); i++) {
        host_constants[totalConstants] = constantsList[i]; totalConstants++;
    }

    //stick placeholders into our normalisation array
    host_normalisations[totalNormalisations] = 1;
    normalIdx = totalNormalisations; totalNormalisations++; 
    host_normalisations[totalNormalisations] = 0;
    totalNormalisations++;

    if(totalParameters >= maxParams)
        throw GooFit::GeneralError("{}: Set too many parameters, GooFit array more than {}. Increase max at compile time with -DGOOFIT_MAXPAR=N.", getName(), maxParams);

    //we rely on GooPdf::set to copy these values, this copyies every PDF which is unnecessary.
    //MEMCPY_TO_SYMBOL(paramIndices, host_indices, totalParams*sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::recursiveSetIndices () {
    //This function always needs to be overloaded.  The key component missing is assigning the function pointer, which is only available
    //in each source file.  Otherwise, we will set the values as follows:

    //This is a helper function if the routine does nothing special.
    //populateArrays ();   
}

__host__ void PdfBase::updateParameters ()
{
    //GOOFIT_TRACE("Update parameters for {}", getName ());

    for (int i = 0; i < parametersList.size (); i++)
        host_parameters[parametersIdx + i + 1] = parametersList[i]->getValue () - parametersList[i]->blind;

    for (int i = 0; i < components.size (); i++)
        components[i]->updateParameters ();

    //we need to memcpy to device.
    MEMCPY(d_parameters, host_parameters, totalParameters*sizeof(fptype), cudaMemcpyHostToDevice);
}

__host__ void PdfBase::populateArrays () {
    //populate all the arrays 
    //GOOFIT_DEBUG("Populating Arrays for {}", getName());
    //GOOFIT_TRACE("host_parameters[{}] = {}", totalParameters, parametersList.size());
    host_parameters[totalParameters] = parametersList.size ();
    parametersIdx = totalParameters; totalParameters++;
    for (int i = 0; i < parametersList.size (); i++)
    {
        //GOOFIT_TRACE("host_parameters[{}] = {}", totalParameters, parametersList[i]->getValue());
        host_parameters[totalParameters] = parametersList[i]->getValue (); totalParameters++;
    }

    //GOOFIT_TRACE("host_constants[{}] = {}", totalConstants, constantsList.size());
    host_constants[totalConstants] = constantsList.size ();
    constantsIdx = totalConstants; totalConstants++;
    for (int i = 0; i < constantsList.size (); i++)
    {
        //GOOFIT_TRACE("host_constants[{}] = {}", totalConstants, constantsList[i]);
        host_constants[totalConstants] = constantsList[i]; totalConstants++;
    }

    //GOOFIT_TRACE("host_observables[{}] = {}", totalObservables, observablesList.size());
    host_observables[totalObservables] = observablesList.size ();
    observablesIdx = totalObservables; totalObservables++;
    for (int i = 0; i < observablesList.size (); i++)
    {
        //GOOFIT_TRACE("host_observables[{}] = {}", totalObservables, observablesList[i]->getValue ());
        host_observables[totalObservables] = observablesList[i]->getValue (); totalObservables++;
    }

    //GOOFIT_TRACE("host_normalisations[{}] = {}", totalNormalisations, 1);
    host_normalisations[totalNormalisations] = 1;
    normalIdx = totalNormalisations++;
    //GOOFIT_TRACE("host_normalisations[{}] = {}", totalNormalisations, 0);
    host_normalisations[totalNormalisations] = 0; totalNormalisations++;

    for (unsigned int i = 0; i < components.size (); i++)
        components[i]->recursiveSetIndices ();

    generateNormRange ();
}

__host__ void PdfBase::setData(std::vector<std::map<Variable*, fptype>>& data) {
    // Old method retained for backwards compatibility

    if(dev_event_array) {
        gooFree(dev_event_array);
        dev_event_array = 0;
    }

    setIndices();
    int dimensions = observablesList.size();
    numEntries = data.size();
    numEvents = numEntries;

    fptype* host_array = new fptype[data.size()*dimensions];

    for(unsigned int i = 0; i < data.size(); ++i) {
        for(Variable*  v : observablesList) {
            if(data[i].find(v) == data[i].end())
                throw GooFit::GeneralError("Variable {} not found", v->getName());
            host_array[i*dimensions + v->getIndex()] = data[i][v];
        }
    }

    gooMalloc((void**) &dev_event_array, dimensions*numEntries*sizeof(fptype));
    MEMCPY(dev_event_array, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
    //MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    delete[] host_array;
}

__host__ void PdfBase::setIndices() {
    int counter = 0;

    //This is unnecessary?
    for(Variable* v : observablesList) {
        v->setIndex(counter++);
    }

    //we should get the same amount after we flatten the tree!
    int checkParams = totalParameters;

    //This is the key, we flatten the tree by re-running through the whole PDF with everything zero'd.
    totalParameters = 0;
    totalConstants = 0;
    totalObservables = 0;
    totalNormalisations = 0;
    num_device_functions = 0;

    //set all associated functions parameters, constants, etc.
    recursiveSetIndices();

    if (checkParams != totalParameters)
        GOOFIT_DEBUG("Error!  checkParams({}) != totalParameters({})", checkParams, totalParameters);
}

__host__ void PdfBase::setData(DataSet* data) {
    if(dev_event_array) {
        gooFree(dev_event_array);
        cudaDeviceSynchronize();
        dev_event_array = 0;
        m_iEventsPerTask = 0;
    }

    setIndices();
    
    UnbinnedDataSet* unbinned_data;
    BinnedDataSet* binned_data;
    
    if((unbinned_data = dynamic_cast<UnbinnedDataSet*>(data))) {
        numEntries = data->getNumEvents();
        numEvents = numEntries;
        
        int dimensions = observablesList.size();

    #ifdef GOOFIT_MPI
        //This fetches our rank and the total number of processes in the MPI call
        int myId, numProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);

        int perTask = numEvents/numProcs;

        //This will track for a given rank where they will start and how far they will go
        int* counts = new int[numProcs];
        int* displacements = new int[numProcs];

        for(int i = 0; i < numProcs - 1; i++)
            counts[i] = perTask;

        counts[numProcs - 1] = numEntries - perTask*(numProcs - 1);

        displacements[0] = 0;

        for(int i = 1; i < numProcs; i++)
            displacements[i] = displacements[i - 1] + counts[i - 1];

    #endif

        fptype* host_array = new fptype[numEntries*dimensions];

    #ifdef GOOFIT_MPI
        //This is an array to track if we need to re-index the observable
        int fixme[observablesList.size()];
        memset(fixme, 0, sizeof(int)*observablesList.size());

        for(int i = 0; i < observables.size(); i++) {
            //We are casting the observable to a CountVariable
            CountingVariable* c = dynamic_cast <CountingVariable*>(observablesList[i]);

            //if it is true re-index
            if(c)
                fixme[i] = 1;
        }

    #endif

        //Transfer into our whole buffer
        for(int i = 0; i < numEntries; ++i) {
            for(Variable* v : observablesList) {
                fptype currVal = unbinned_data->getValue(v, i);
                host_array[i*dimensions + v->getIndex()] = currVal;
            }
        }

    #ifdef GOOFIT_MPI

        //We will go through all of the events and re-index if appropriate
        for(int i = 1; i < numProcs; i++) {
            for(int j = 0; j < counts[i]; j++) {
                for(int k = 0; k < dimensions; k++) {
                    if(fixme[k] > 0)
                        host_array[(j + displacements[i])*dimensions + k] = float (j);
                }
            }
        }

        int mystart = displacements[myId];
        int myend = mystart + counts[myId];
        int mycount = myend - mystart;

        gooMalloc((void**) &dev_event_array, dimensions*mycount*sizeof(fptype));
        MEMCPY(dev_event_array, host_array + mystart*dimensions, dimensions*mycount*sizeof(fptype), cudaMemcpyHostToDevice);
        //MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
        delete[] host_array;

        setNumPerTask(this, mycount);

        delete []counts;
        delete []displacements;
    #else
        gooMalloc((void**) &dev_event_array, dimensions*numEntries*sizeof(fptype));
        MEMCPY(dev_event_array, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
        //MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
        delete[] host_array;
    #endif
    } else if ((binned_data = dynamic_cast<BinnedDataSet*>(data))) {

     
        numEvents = 0;
        numEntries = binned_data->getNumBins();
        int dimensions = 2 + observablesList.size(); // Bin center (x,y, ...), bin value, and bin volume.

        if(!fitControl->binnedFit())
            setFitControl(new BinnedNllFit());

    #ifdef GOOFIT_MPI
        //This fetches our rank and the total number of processes in the MPI call
        int myId, numProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);

        int perTask = numEvents/numProcs;

        //This will track for a given rank where they will start and how far they will go
        int* counts = new int[numProcs];
        int* displacements = new int[numProcs];

        for(int i = 0; i < numProcs - 1; i++)
            counts[i] = perTask;

        counts[numProcs - 1] = numEntries - perTask*(numProcs - 1);

        displacements[0] = 0;

        for(int i = 1; i < numProcs; i++)
            displacements[i] = displacements[i - 1] + counts[i - 1];

    #endif

        fptype* host_array = new fptype[numEntries*dimensions];

    #ifdef GOOFIT_MPI
        //This is an array to track if we need to re-index the observable
        int fixme[observables.size()];
        memset(fixme, 0, sizeof(int)*observablesList.size());

        for(int i = 0; i < observablesList.size(); i++) {
            //We are casting the observable to a CountVariable
            CountingVariable* c = dynamic_cast <CountingVariable*>(observablesList[i]);

            //if it is true re-index
            if(c)
                fixme[i] = 1;
        }

    #endif

        for(unsigned int i = 0; i < numEntries; ++i) {
            for(Variable* v : observablesList) {
                host_array[i*dimensions + v->getIndex()] = binned_data->getBinCenter(v, i);
            }

            host_array[i*dimensions + observablesList.size() + 0] = binned_data->getBinContent(i);
            host_array[i*dimensions + observablesList.size() + 1] = fitControl->binErrors() ? binned_data->getBinError(i) : binned_data->getBinVolume(i);
            numEvents += binned_data->getBinContent(i);
        }

    #ifdef GOOFIT_MPI

        //We will go through all of the events and re-index if appropriate
        for(int i = 1; i < numProcs; i++) {
            for(int j = 0; j < counts[j]; j++) {
                for(int k = 0; k < dimensions; k++) {
                    if(fixme[k] > 0)
                        host_array[(j + displacements[i])*dimensions + k] = float (j);
                }
            }
        }

        int mystart = displacements[myId];
        int myend = mystart + counts[myId];
        int mycount = myend - mystart;

        gooMalloc((void**) &dev_event_array, dimensions*mycount*sizeof(fptype));
        MEMCPY(dev_event_array, host_array + mystart*dimensions, dimensions*mycount*sizeof(fptype), cudaMemcpyHostToDevice);
        //MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
        delete[] host_array;

        setNumPerTask(this, mycount);

        delete []counts;
        delete []displacements;
    #else
        gooMalloc((void**) &dev_event_array, dimensions*numEntries*sizeof(fptype));
        MEMCPY(dev_event_array, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
        //MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
        delete[] host_array;
    #endif
    } else
        throw GooFit::GeneralError("Dataset must be binned or unbinned!");
}

__host__ void PdfBase::generateNormRange() {
    if(normRanges)
        gooFree(normRanges);

    gooMalloc((void**) &normRanges, 3*observablesList.size()*sizeof(fptype));

    fptype* host_norms = new fptype[3*observablesList.size()];
    int counter = 0; // Don't use index in this case to allow for, eg,

    // a single observable whose index is 1; or two observables with indices
    // 0 and 2. Make one array per functor, as opposed to variable, to make
    // it easy to pass MetricTaker a range without worrying about which parts
    // to use.
    for(Variable* v : observablesList) {
        host_norms[3*counter+0] = v->getLowerLimit();
        host_norms[3*counter+1] = v->getUpperLimit();
        host_norms[3*counter+2] = integrationBins > 0 ? integrationBins : v->getNumBins();
        GOOFIT_TRACE("host_norms[{}] = {}", 3*counter + 0, host_norms[3*counter + 0]);
        GOOFIT_TRACE("host_norms[{}] = {}", 3*counter + 1, host_norms[3*counter + 1]);
        GOOFIT_TRACE("host_norms[{}] = {}", 3*counter + 2, host_norms[3*counter + 2]);
        counter++;
    }

    MEMCPY(normRanges, host_norms, 3*observablesList.size()*sizeof(fptype), cudaMemcpyHostToDevice);
    delete[] host_norms;
}

void PdfBase::clearCurrentFit() {
    totalParameters = 0;
    gooFree(dev_event_array);
    dev_event_array = 0;
}

__host__ void PdfBase::printProfileInfo(bool topLevel) {
#ifdef PROFILING

    if(topLevel) {
        cudaError_t err = MEMCPY_FROM_SYMBOL(host_timeHist, timeHistogram, 10000*sizeof(fptype), 0);

        if(cudaSuccess != err) {
            std::cout << "Error on copying timeHistogram: " << cudaGetErrorString(err) << std::endl;
            return;
        }

        std::cout << getName() << " : " << getFunctionIndex() << " " << host_timeHist[100*getFunctionIndex() +
                                         getParameterIndex()] << std::endl;

        for(unsigned int i = 0; i < components.size(); ++i) {
            components[i]->printProfileInfo(false);
        }
    }

#endif
}



gooError gooMalloc(void** target, size_t bytes) {
#if THRUST_DEVICE_SYSTEM!=THRUST_DEVICE_SYSTEM_CUDA
    target[0] = malloc(bytes);

    if(target[0])
        return gooSuccess;
    else
        return gooErrorMemoryAllocation;

#else
    return (gooError) cudaMalloc(target, bytes);
#endif
}

gooError gooFree(void* ptr) {
#if THRUST_DEVICE_SYSTEM!=THRUST_DEVICE_SYSTEM_CUDA
    free(ptr);
    return gooSuccess;
#else
    return (gooError) cudaFree(ptr);
#endif
}
} // namespace GooFit

