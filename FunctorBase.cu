#include "FunctorBase.hh"

// This is code that belongs to the FunctorBase class, that is, 
// it is common across all implementations. But it calls on device-side
// functions, and due to the nvcc translation-unit limitations, it cannot
// sit in its own object file; it must go in the CUDAglob.cu. So it's
// off on its own in this inline-cuda file, which ThrustPdfFunctor.cu 
// should include. 

#ifdef CUDAPRINT
__host__ void FunctorBase::copyParams (const std::vector<double>& pars) const {
  if (host_callnumber < 1) {
    std::cout << "Copying parameters: " << (long long) cudaArray << " ";
  }
  for (unsigned int i = 0; i < pars.size(); ++i) {
    host_params[i] = pars[i]; 
    
    if (host_callnumber < 1) {
      std::cout << pars[i] << " ";
    }
    
    if (isnan(host_params[i])) {
      std::cout << " agh, NaN, die " << i << std::endl;
      abortWithCudaPrintFlush(__FILE__, __LINE__, "NaN in parameter"); 
    }
  }
  
  if (host_callnumber < 1) {
    std::cout << std::endl; 
  }
  cudaMemcpyToSymbol(cudaArray, host_params, pars.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
}
#else 
__host__ void FunctorBase::copyParams (const std::vector<double>& pars) const {
  for (unsigned int i = 0; i < pars.size(); ++i) {
    host_params[i] = pars[i]; 
    
    if (isnan(host_params[i])) {
      std::cout << " agh, parameter is NaN, die " << i << std::endl;
      abortWithCudaPrintFlush(__FILE__, __LINE__, "NaN in parameter"); 
    }
  }

  cudaMemcpyToSymbol(cudaArray, host_params, pars.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
}
#endif

__host__ void FunctorBase::copyParams () {
  // Copies values of Variable objects
  std::set<Variable*> pars;
  getParameters(pars); 
  std::vector<double> values; 
  //if (cpuDebug & 1) std::cout << "Copying parameters: "; 
  for (std::set<Variable*>::iterator v = pars.begin(); v != pars.end(); ++v) {
    int index = (*v)->getIndex(); 
    if (index >= (int) values.size()) values.resize(index + 1);
    values[index] = (*v)->value;
    //if (cpuDebug & 1) std::cout << "(" << (*v)->index << " " << (*v)->value << ") ";
    //if (isnan((*v)->value)) std::cout << "Major issue, " << (*v)->name << " is NaN.\n"; 
  }
  //if (cpuDebug & 1) std::cout << std::endl;
  copyParams(values); 
}

__host__ void FunctorBase::copyNormFactors () const {
  cudaMemcpyToSymbol(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  cudaDeviceSynchronize(); // Ensure normalisation integrals are finished
}

__host__ void FunctorBase::initialiseIndices (std::vector<unsigned int> pindices) {
  // Structure of the individual index array: Number of parameters, then the indices
  // requested by the subclass (which will be interpreted by the subclass kernel), 
  // then the number of observables, then the observable indices. Notice that the
  // observable indices are not set until 'setIndices' is called, usually from setData;
  // here we only reserve space for them by setting totalParams. 
  // This is to allow index sharing between PDFs - all the PDFs must be constructed 
  // before we know what observables exist. 

  if (totalParams + pindices.size() >= maxParams) {
    std::cout << "Major problem with pindices size: " << totalParams << " + " << pindices.size() << " >= " << maxParams << std::endl; 
  }

  assert(totalParams + pindices.size() < maxParams); 
  host_indices[totalParams] = pindices.size(); 
  for (int i = 1; i <= host_indices[totalParams]; ++i) {
    host_indices[totalParams+i] = pindices[i-1]; 
  }
  host_indices[totalParams + pindices.size() + 1] = observables.size(); 
  
  parameters = totalParams;
  totalParams += (2 + pindices.size() + observables.size()); 
  /*
  std::cout << "host_indices after " << getName() << " initialisation : ";
  for (int i = 0; i < totalParams; ++i) {
    std::cout << host_indices[i] << " ";
  }
  
  std::cout << " | " 
	    << parameters << " " 
	    << totalParams << " " 
	    << cudaArray << " " 
	    << paramIndices << " "
	    << std::endl; 
  */
  cudaMemcpyToSymbol(paramIndices, host_indices, totalParams*sizeof(unsigned int), 0, cudaMemcpyHostToDevice); 
}

__host__ void FunctorBase::setData (std::vector<std::map<Variable*, fptype> >& data) {
  // Old method retained for backwards compatibility 

  if (cudaDataArray) {
    cudaFree(cudaDataArray);
    cudaDataArray = 0; 
  }

  setIndices();
  int dimensions = observables.size();
  numEntries = data.size();
  numEvents = numEntries; 
  
  fptype* host_array = new fptype[data.size()*dimensions];
  for (unsigned int i = 0; i < data.size(); ++i) {
    for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
      assert(data[i].find(*v) != data[i].end()); 
      host_array[i*dimensions + (*v)->index] = data[i][*v]; 
    }
  }

  cudaMalloc((void**) &cudaDataArray, dimensions*numEntries*sizeof(fptype)); 
  cudaMemcpy(cudaDataArray, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  delete[] host_array; 
}

__host__ void FunctorBase::recursiveSetIndices () {
  for (unsigned int i = 0; i < components.size(); ++i) {
    components[i]->recursiveSetIndices(); 
  }

  int numParams = host_indices[parameters]; 
  int counter = 0; 
  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    host_indices[parameters + 2 + numParams + counter] = (*v)->index; 
    //std::cout << getName() << " set index of " << (*v)->name << " to " << (*v)->index << " " << (parameters + 2 + numParams + counter) << std::endl; 
    counter++; 
  }  
  generateNormRange(); 
}

__host__ void FunctorBase::setIndices () {
  int counter = 0; 
  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    (*v)->index = counter++; 
  }
  recursiveSetIndices(); 
  cudaMemcpyToSymbol(paramIndices, host_indices, totalParams*sizeof(unsigned int), 0, cudaMemcpyHostToDevice); 

  //std::cout << "host_indices after " << getName() << " observable setIndices : ";
  //for (int i = 0; i < totalParams; ++i) {
  //std::cout << host_indices[i] << " ";
  //}
  //std::cout << std::endl; 

}

__host__ void FunctorBase::setData (UnbinnedDataSet* data) {
  if (cudaDataArray) {
    cudaFree(cudaDataArray);
    cudaDeviceSynchronize();
    cudaDataArray = 0; 
  }

  setIndices();
  int dimensions = observables.size();
  numEntries = data->getNumEvents(); 
  numEvents = numEntries; 
  
  fptype* host_array = new fptype[numEntries*dimensions];
  for (int i = 0; i < numEntries; ++i) {
    for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
      fptype currVal = data->getValue((*v), i);
      host_array[i*dimensions + (*v)->index] = currVal; 
      //if (60 == i) std::cout << "Setting data point " << (*v)->name << " " << currVal << " " << i*dimensions + (*v)->index << " " << (*v)->index << std::endl; 
    }
  }

  //std::cout << "Trying to allocate " << (numEntries*dimensions*sizeof(fptype)) << " = " << numEntries << " * " << dimensions << " * " << sizeof(fptype) << " bytes.\n"; 
  cudaMalloc((void**) &cudaDataArray, dimensions*numEntries*sizeof(fptype)); 
  //cudaDeviceSynchronize();
  //std::cout << "Allocated to " << cudaDataArray << std::endl; 
  cudaMemcpy(cudaDataArray, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  delete[] host_array; 
}

__host__ void FunctorBase::setData (BinnedDataSet* data) { 
  if (cudaDataArray) { 
    cudaFree(cudaDataArray);
    cudaDataArray = 0; 
  }

  setIndices();
  numEvents = 0; 
  numEntries = data->getNumBins(); 
  int dimensions = 2 + observables.size(); // Bin center (x,y, ...), bin value, and bin volume. 
  if (!fitControl->binnedFit()) setFitControl(new BinnedNllFit()); 

  fptype* host_array = new fptype[numEntries*dimensions]; 

  for (unsigned int i = 0; i < numEntries; ++i) {
    for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
      host_array[i*dimensions + (*v)->index] = data->getBinCenter((*v), i); 
    }

    host_array[i*dimensions + observables.size() + 0] = data->getBinContent(i);
    host_array[i*dimensions + observables.size() + 1] = fitControl->binErrors() ? data->getBinError(i) : data->getBinVolume(i); 
    numEvents += data->getBinContent(i);
  }

  cudaMalloc(&cudaDataArray, dimensions*numEntries*sizeof(fptype)); 
  cudaMemcpy(cudaDataArray, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  delete[] host_array; 
}

/*
__host__ void FunctorBase::setBinnedData (std::vector<std::pair<std::map<Variable*, fptype>, fptype> >& data) {
  if (cudaDataArray) { 
    cudaFree(cudaDataArray);
    cudaDataArray = 0; 
  }

  numEvents = 0; 
  numEntries = data.size(); 
  int dimensions = 2 + observables.size(); // Bin center (x,y, ...), bin value, and bin volume. 

  fptype* host_array = new fptype[numEntries*dimensions]; 
  int count = 0; 


  for (unsigned int i = 0; i < numEntries; ++i) {
    fptype binVolume = 1; 
    for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
      binVolume *= ((*v)->upperlimit - (*v)->lowerlimit) / (*v)->numbins; 
      //      assert(data[i].first.find(*v) != data[i].first.end()); 
      host_array[count++] = data[i].first[*v]; 
      //host_array[count++] = (*v)->lowerlimit + (i + 0.5) * (((*v)->upperlimit - (*v)->lowerlimit) / (*v)->numbins);
    }

    host_array[count++] = data[i].second; 
    host_array[count++] = binVolume; 
    numEvents += data[i].second; 
  }

  cudaMalloc(&cudaDataArray, dimensions*numEntries*sizeof(fptype)); 
  cudaMemcpy(cudaDataArray, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  delete[] host_array; 

  std::cout << "Line " << __LINE__ << "\n";
  
  if (!fitControl->binnedFit()) setFitControl(new BinnedNllFit()); 
}
*/

__host__ void FunctorBase::generateNormRange () {
  if (normRanges) cudaFree(normRanges);
  cudaMalloc((void**) &normRanges, 3*observables.size()*sizeof(fptype));
  
  fptype* host_norms = new fptype[3*observables.size()];
  int counter = 0; // Don't use index in this case to allow for, eg, 
  // a single observable whose index is 1; or two observables with indices
  // 0 and 2. Make one array per functor, as opposed to variable, to make
  // it easy to pass MetricTaker a range without worrying about which parts
  // to use. 
  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    host_norms[3*counter+0] = (*v)->lowerlimit;
    host_norms[3*counter+1] = (*v)->upperlimit;
    host_norms[3*counter+2] = integrationBins > 0 ? integrationBins : (*v)->numbins;
    counter++; 

    //std::cout << "Norm range " << getName() << ": " << (*v)->name << " " << (*v)->index << " " << host_norms[3*(*v)->index+0] << " " << host_norms[3*(*v)->index+1] << " " << host_norms[3*(*v)->index+2] << std::endl; 
  }

  cudaMemcpy(normRanges, host_norms, 3*observables.size()*sizeof(fptype), cudaMemcpyHostToDevice);
  //std::cout << "Done with copy\n";
  delete[] host_norms; 
  //std::cout << "Done with delete\n"; 
}

void FunctorBase::clearCurrentFit () {
  totalParams = 0; 
  cudaFree(cudaDataArray);
  cudaDataArray = 0; 
}

__host__ void FunctorBase::printProfileInfo (bool topLevel) {
#ifdef PROFILING
  if (topLevel) {
    cudaError_t err = cudaMemcpyFromSymbol(host_timeHist, timeHistogram, 10000*sizeof(fptype), 0);
    if (cudaSuccess != err) {
      std::cout << "Error on copying timeHistogram: " << cudaGetErrorString(err) << std::endl;
      return;
    }
    /*
    std::cout << "Time histogram: "; 
    for (int i = 0; i < num_device_functions; ++i) {
      std::cout << host_timeHist[i] << " ";
    }
    std::cout << std::endl;
    */ 
  }

  std::cout << getName() << " : " << getFunctionIndex() << " " << host_timeHist[100*getFunctionIndex() + getParameterIndex()] << std::endl; 
  for (unsigned int i = 0; i < components.size(); ++i) {
    components[i]->printProfileInfo(false); 
  }
#endif
}
