#include "goofit/PDFs/combine/ConvolutionPdf.h"
#include "goofit/Error.h"
#include "goofit/Variable.h"

#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>

namespace GooFit {

int totalConvolutions = 0;

#define CONVOLUTION_CACHE_SIZE 512
// I would really like this to be 1024, which is the maximum number of threads
// in a block (for compute capability 2.0 and up). Unfortunately this causes
// the program to hang, presumably because there isn't enough memory and something
// goes wrong. So... 512 should be enough for anyone, right?

// Need multiple working spaces for the case of several convolutions in one PDF.
__constant__ fptype *dev_modWorkSpace[100];
__constant__ fptype *dev_resWorkSpace[100];

// Number which transforms model range (x1, x2) into resolution range (x1 - maxX, x2 - minX).
// It is equal to the maximum possible value of x0, ie maxX, in bins.
__constant__ int modelOffset[100];

__device__ fptype device_ConvolvePdfs(fptype* evt, ParameterContainer &pc) {
    int id = pc.constants[pc.constantIdx + 2];

    fptype ret     = 0;
    fptype loBound = pc.constants[pc.constantIdx + 3];//RO_CACHE(pc.constants[pc.constantIdx + 2]);
    fptype hiBound = pc.constants[pc.constantIdx + 4];//RO_CACHE(pc.constants[pc.constantIdx + 3]);
    fptype step    = pc.constants[pc.constantIdx + 5];//RO_CACHE(pc.constants[pc.constantIdx + 4]);
    fptype x0      = evt[id];
    int workSpaceIndex = pc.constants[pc.constantIdx + 6];

    auto numbins = static_cast<int>(floor((hiBound - loBound) / step + 0.5));

    fptype lowerBoundOffset = loBound / step;
    lowerBoundOffset -= floor(lowerBoundOffset);
    auto offsetInBins = static_cast<int>(floor(x0 / step - lowerBoundOffset));

    // Brute-force calculate integral M(x) * R(x - x0) dx
    int offset = RO_CACHE(modelOffset[workSpaceIndex]);

    for(int i = 0; i < numbins; ++i) {
        fptype model = RO_CACHE(dev_modWorkSpace[workSpaceIndex][i]);
        fptype resol = RO_CACHE(dev_resWorkSpace[workSpaceIndex][i + offset - offsetInBins]);

        ret += model * resol;
    }

    //What the heck is this?
    //When the convolution is flattened, the model PDF and resolution PDF are maintained and linked.
    //These increments will skip the model & resolution PDF.  We will need to determine if we need a
    //skip all chilren also event.  to be determined...
    pc.incrementIndex (1, 0, 6, 0, 1);

    fptype norm1 = pc.normalisations[pc.normalIdx + 1];
    ret *= norm1;

    pc.incrementIndex ();

    fptype norm2 = pc.normalisations[pc.normalIdx + 1];
    ret *= norm2;

    pc.incrementIndex ();

    return ret;
}

__device__ fptype device_ConvolveSharedPdfs(fptype* evt, ParameterContainer &pc) {
    fptype ret     = 0;
    fptype loBound = pc.constants[pc.constantIdx + 3];
    fptype hiBound = pc.constants[pc.constantIdx + 4];
    fptype step    = pc.constants[pc.constantIdx + 5];
    fptype x0      = evt[0];
    unsigned int workSpaceIndex = pc.constants[pc.constantIdx + 1];
    unsigned int numOthers      = pc.constants[pc.constantIdx + 2] + 1; // +1 for this PDF.

    auto numbins = static_cast<int>(floor((hiBound - loBound) / step + 0.5));

    fptype lowerBoundOffset = loBound / step;
    lowerBoundOffset -= floor(lowerBoundOffset);
    auto offsetInBins = static_cast<int>(floor(x0 / step - lowerBoundOffset));

    // Brute-force calculate integral M(x) * R(x - x0) dx
    __shared__ fptype modelCache[CONVOLUTION_CACHE_SIZE];
    // Don't try to shared-load more items than we have threads.
    int numToLoad
        = thrust::minimum<unsigned int>()(CONVOLUTION_CACHE_SIZE / numOthers, static_cast<unsigned int>(BLOCKDIM));

    for(int i = 0; i < numbins; i += numToLoad) {
        // This code avoids this problem: If event 0 is in workspace 0, and
        // event 1 is in workspace 1, then threads 0 and 1 cannot fill the
        // same modelCache, because they are drawing from different models!
        // What's more, it's not guaranteed that the whole modelCache will
        // be filled, because the threads can be all over the place in terms
        // of their workspace.

        if(THREADIDX < numToLoad) {
            for(unsigned int w = 0; w < numOthers; ++w) {
                unsigned int wIndex = pc.constants[pc.constantIdx + 5 + w];
                modelCache[w * numToLoad + THREADIDX]
                    = (i + THREADIDX < numbins) ? dev_modWorkSpace[wIndex][i + THREADIDX] : 0;
            }
        }

        // This is slightly dangerous. If you have, for example,
        // a MappedPdf in which one of the other functions
        // is not a convolution, *you're gonna have a bad time*.
        // It's unfortunately up to the user to ensure that every
        // event will hit this point, by making sure that if *any*
        // event has a convolution, *all* events do.
        THREAD_SYNCH

        for(int j = 0; j < numToLoad; ++j) {
            if(i + j >= numbins)
                break;

            fptype model = modelCache[workSpaceIndex * numToLoad + j];
            fptype resol = (model != 0)
                               ? dev_resWorkSpace[workSpaceIndex][i + j + modelOffset[workSpaceIndex] - offsetInBins]
                               : 0;
            ret += model * resol;
        }

        // If we don't sync at this point, some of the warps can *run ahead*
        // and fill up their parts of modelCache with the *next* set of numbers,
        // which means that the warps still working on the sum get the wrong
        // numbers!
        THREAD_SYNCH
    }

    ret *= pc.normalisations[pc.normalIdx + 1];
    ret *= pc.normalisations[pc.normalIdx + 2];

    return ret;
}

__device__ device_function_ptr ptr_to_ConvolvePdfs       = device_ConvolvePdfs;
__device__ device_function_ptr ptr_to_ConvolveSharedPdfs = device_ConvolveSharedPdfs;

ConvolutionPdf::ConvolutionPdf(std::string n, Variable *x, GooPdf *m, GooPdf *r)
    : GooPdf(x, n)
    , model(m)
    , resolution(r)
    , host_iConsts(nullptr)
    , modelWorkSpace(nullptr)
    , resolWorkSpace(nullptr)
    , workSpaceIndex(0) {
    // Constructor for convolution without cooperative
    // loading of model cache. This is slow, but conceptually
    // simple.

    observablesList = getObservables();

    //reserving index for x
    constantsList.push_back(observablesList.size());
    constantsList.push_back (0);

    components.push_back(model);
    components.push_back(resolution);

    // Indices stores (function index)(parameter index) doublet for model and resolution function.
    std::vector<unsigned int> paramIndices;
    //paramIndices.push_back(model->getFunctionIndex());
    //paramIndices.push_back(model->getParameterIndex());
    //paramIndices.push_back(resolution->getFunctionIndex());
    //paramIndices.push_back(resolution->getParameterIndex());
    //paramIndices.push_back(registerConstants(3));
    //paramIndices.push_back(workSpaceIndex = totalConvolutions++);

    constantsList.push_back (-10);
    constantsList.push_back (10);
    constantsList.push_back (0.01);
    constantsList.push_back (workSpaceIndex = totalConvolutions++);

    GET_FUNCTION_ADDR(ptr_to_ConvolvePdfs);
    initialize(paramIndices);
    setIntegrationConstants(-10, 10, 0.01);

    ConvolveType = 0;
}

ConvolutionPdf::ConvolutionPdf(std::string n, Variable *x, GooPdf *m, GooPdf *r, unsigned int numOthers)
    : GooPdf(x, n)
    , model(m)
    , resolution(r)
    , host_iConsts(nullptr)
    , modelWorkSpace(nullptr)
    , resolWorkSpace(nullptr)
    , workSpaceIndex(0) {
    // Constructor for the case when several convolutions
    // can run in parallel; for example, a map where the targets
    // are convolutions. (Make sure that *all* the targets
    // are convolutions! Otherwise it will hang on the syncthreads
    // call.) In this case the cooperative loading needs
    // to initialize all the shared memory workspaces, and needs to know
    // how many such workspaces there are, and which global workspaces
    // to draw on. NB! To use cooperative loading in the case of a
    // single function, just use numOthers = 0.

    //reserving index for x
    constantsList.push_back (0);

    components.push_back(model);
    components.push_back(resolution);

    // Indices stores (function index)(parameter index) doublet for model and resolution function.
    std::vector<unsigned int> paramIndices;
    //paramIndices.push_back(model->getFunctionIndex());
    //paramIndices.push_back(model->getParameterIndex());
    //paramIndices.push_back(resolution->getFunctionIndex());
    //paramIndices.push_back(resolution->getParameterIndex());
    //paramIndices.push_back(registerConstants(3));
    //paramIndices.push_back(workSpaceIndex = totalConvolutions++);
    //paramIndices.push_back(numOthers);
    constantsList.push_back (0);

    constantsList.push_back (-10);
    constantsList.push_back (10);
    constantsList.push_back (0.01);
    constantsList.push_back (workSpaceIndex = totalConvolutions++);

    if(0 == numOthers)
        constantsList.push_back(workSpaceIndex);
    else {
        properlyInitialised = false;

        for(unsigned int i = 0; i < numOthers + 1; ++i) { // Notice extra space for this PDF's index.
            // Fill in later - must be done before setData call.
            constantsList.push_back(0);
        }
    }

    if(numOthers > CONVOLUTION_CACHE_SIZE)
        throw GooFit::GeneralError(
            "numOthers {} must be not be more than the cache size {}", numOthers, CONVOLUTION_CACHE_SIZE);

    GET_FUNCTION_ADDR(ptr_to_ConvolveSharedPdfs);
    initialize(paramIndices);
    setIntegrationConstants(-10, 10, 0.01);

    ConvolveType = 1;
}

__host__ void ConvolutionPdf::recursiveSetIndices () {
    if (ConvolveType == 0) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_ConvolvePdfs");
        GET_FUNCTION_ADDR(ptr_to_ConvolvePdfs);
    }
    else if (ConvolveType == 1) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_ConvolveSharedPdfs");
        GET_FUNCTION_ADDR(ptr_to_ConvolveSharedPdfs);
    }

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions ++;

    populateArrays ();
}

__host__ void ConvolutionPdf::setIntegrationConstants(fptype lo, fptype hi, fptype step) {
    if(!host_iConsts) {
        host_iConsts = new fptype[6];
        gooMalloc(reinterpret_cast<void **>(&dev_iConsts), 6 * sizeof(fptype));
    }


    host_iConsts[0] = lo;
    host_iConsts[1] = hi;
    host_iConsts[2] = step;

    constantsList[2] = lo;
    constantsList[3] = hi;
    constantsList[4] = step;
    //MEMCPY_TO_SYMBOL(functorConstants, host_iConsts, 3 * sizeof(fptype), cIndex * sizeof(fptype), cudaMemcpyHostToDevice);

    if(modelWorkSpace) {
        delete modelWorkSpace;
        delete resolWorkSpace;
    }

    auto numbins = static_cast<int>(floor((host_iConsts[1] - host_iConsts[0]) / step + 0.5));
    // Different format for integration range!
    modelWorkSpace = new thrust::device_vector<fptype>(numbins);

    // We will do integral from x1 to x2 of M(x)*R(x - x0) dx.
    // So we need to cache the values of M from x1 to x2, which is given
    // by the integration range. But R must be cached from x1-maxX to
    // x2-minX, and the min and max are given by the dependent variable.
    // However, the step must be the same as for the model, or the binning
    // will get out of sync.
    Variable *dependent = *(observablesList.begin());

    host_iConsts[2] = numbins;
    host_iConsts[3] = (host_iConsts[0] - dependent->getUpperLimit());
    host_iConsts[4] = (host_iConsts[1] - dependent->getLowerLimit());

    numbins         = static_cast<int>(floor((host_iConsts[4] - host_iConsts[3]) / step + 0.5));
    host_iConsts[5] = numbins;
    MEMCPY(dev_iConsts, host_iConsts, 6 * sizeof(fptype), cudaMemcpyHostToDevice);
    //MEMCPY(convolutionConstants, host_iConsts, 6 * sizeof(fptype), cudaMemcpyHostToDevice);
    resolWorkSpace = new thrust::device_vector<fptype>(numbins);

    int offset = dependent->getUpperLimit() / step;
    MEMCPY_TO_SYMBOL(modelOffset, &offset, sizeof(int), workSpaceIndex * sizeof(int), cudaMemcpyHostToDevice);

    fptype *dev_address[1];
    dev_address[0] = (&((*modelWorkSpace)[0])).get();
    MEMCPY_TO_SYMBOL( dev_modWorkSpace, dev_address, sizeof(fptype *), workSpaceIndex * sizeof(fptype *), cudaMemcpyHostToDevice);
    dev_address[0] = (&((*resolWorkSpace)[0])).get();
    MEMCPY_TO_SYMBOL( dev_resWorkSpace, dev_address, sizeof(fptype *), workSpaceIndex * sizeof(fptype *), cudaMemcpyHostToDevice);
}

__host__ void ConvolutionPdf::registerOthers(std::vector<ConvolutionPdf *> others) {
    unsigned int numExpectedOthers = constantsList[constantsIdx + 5] + 1;

    if(numExpectedOthers != others.size()) {
        std::cout << "Problem: " << getName() << " initialized with " << others.size() << " other PDFs, expected "
                  << numExpectedOthers << " (including itself). Returning without initialization.\n";
        return;
    }

    bool foundSelf = false;

    for(unsigned int i = 0; i < others.size(); ++i) {
        ConvolutionPdf *curr = others[i];

        if(curr == this)
            foundSelf = true;

        constantsList[5 + i] = curr->workSpaceIndex;
    }

    if(!foundSelf) {
        std::cout << "Problem: " << getName()
                  << " initialized with list that did not include itself. Returning without initialisation.\n";
        return;
    }

    properlyInitialised = true;
}

__host__ fptype ConvolutionPdf::normalize() const {
    // First set normalisation factors to one so we can evaluate convolution without getting zeroes
    recursiveSetNormalisation(fptype(1.0));

    //we need to update the normal here, as values are used at this point.
    MEMCPY_TO_SYMBOL(d_normalisations, host_normalisations, totalNormalisations*sizeof(fptype), 0, cudaMemcpyHostToDevice);

    // Next recalculate functions at each point, in preparation for convolution integral
    thrust::constant_iterator<fptype *> arrayAddress(dev_iConsts);
    thrust::constant_iterator<int> eventSize(1);
    thrust::counting_iterator<int> binIndex(0);

    if(model->parametersChanged()) {
        // Calculate model function at every point in integration space
        MetricTaker modalor(model, getMetricPointer("ptr_to_Eval"));
        thrust::transform(
            thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress)),
            thrust::make_zip_iterator(thrust::make_tuple(binIndex + modelWorkSpace->size(), eventSize, arrayAddress)),
            modelWorkSpace->begin(),
            modalor);
        cudaDeviceSynchronize();
    }

    if(resolution->parametersChanged()) {
        // Same for resolution function.
        thrust::constant_iterator<fptype *> arrayAddress2(dev_iConsts + 3);
        thrust::constant_iterator<int> resFunc (resolution->getFunctionIndex ());
        MetricTaker resalor(resolution, getMetricPointer("ptr_to_Eval"));
        thrust::transform(
            thrust::make_zip_iterator(thrust::make_tuple(binIndex, eventSize, arrayAddress2)),
            thrust::make_zip_iterator(thrust::make_tuple(binIndex + resolWorkSpace->size(), eventSize, arrayAddress2)),
            resolWorkSpace->begin(),
            resalor);
    }

    // Then return usual integral
    fptype ret = GooPdf::normalize();
    return ret;
}

} // namespace GooFit
