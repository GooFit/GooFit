#include "CompositePdf.hh"

EXEC_TARGET fptype device_Composite (fptype* evt, fptype* p, unsigned int* indices) {
  unsigned int coreFcnIndex  = RO_CACHE(indices[1]);
  unsigned int coreParIndex  = RO_CACHE(indices[2]);
  unsigned int shellFcnIndex = RO_CACHE(indices[3]);
  unsigned int shellParIndex = RO_CACHE(indices[4]);

  // NB, not normalising core function, it is not being used as a PDF. 
  //fptype coreValue = (*(reinterpret_cast<device_function_ptr>(device_function_table[coreFcnIndex])))(evt, cudaArray, paramIndices+coreParIndex);
  fptype coreValue = callFunction(evt, coreFcnIndex, coreParIndex);

  unsigned int* shellParams = paramIndices + shellParIndex; 
  unsigned int numShellPars = shellParams[0];
  unsigned int shellObsIndex = shellParams[2 + numShellPars];

  fptype fakeEvt[10]; // Allow plenty of space in case events are large. 
  fakeEvt[shellObsIndex] = coreValue; 

  // Don't normalise shell either, since we don't know what composite function is being used for. 
  // It may not be a PDF. Normalising at this stage would be presumptuous. 
  //fptype ret = (*(reinterpret_cast<device_function_ptr>(device_function_table[shellFcnIndex])))(fakeEvt, cudaArray, shellParams); 
  fptype ret = callFunction(fakeEvt, shellFcnIndex, shellParIndex); 

  //if (0 == THREADIDX) 
  //printf("Composite: %f %f %f %f %f %f\n", evt[4], evt[5], evt[6], evt[7], coreValue, ret); 

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_Composite = device_Composite; 

__host__ CompositePdf::CompositePdf (std::string n, PdfBase* core, PdfBase* shell) 
  : GooPdf(0, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(core->getFunctionIndex());
  pindices.push_back(core->getParameterIndex());
  pindices.push_back(shell->getFunctionIndex());
  pindices.push_back(shell->getParameterIndex());

  // Add as components so that observables and parameters will be registered. 
  components.push_back(core);
  components.push_back(shell);

  GET_FUNCTION_ADDR(ptr_to_Composite);
  initialise(pindices); 
}

__host__ fptype CompositePdf::normalise () const {
  recursiveSetNormalisation(1.0); 

  // Note: Core is not normalised in composite calculation, 
  // because it is not a PDF, 
  // it is just a plain old function;
  // it can take any value. 
  // Shell needn't be normalised either,
  // because we don't know that the composite
  // will be used as a PDF; and if it is, the
  // normalisation should be applied at the level
  // of whatever calls the composite.
  // However: These functions may appear elsewhere
  // in the full function, and perhaps need to 
  // be normalised there. Consequently, we
  // normalise them even though the information
  // may not be used. 

  for (std::vector<PdfBase*>::const_iterator c = components.begin(); c != components.end(); ++c) {
    (*c)->normalise(); 
  }
  
  // Normalise composite in the usual binned-integral way. 
  return GooPdf::normalise(); 

}
