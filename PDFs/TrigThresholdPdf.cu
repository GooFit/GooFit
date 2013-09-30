#include "TrigThresholdPdf.hh"

EXEC_TARGET fptype threshCalc (fptype distance, fptype linConst) {
  fptype ret = (distance > fptype(0.5) ? fptype(1) : (linConst + (1 - linConst) * SIN(distance * fptype(3.14159265)))); 
  return ret; 
}

EXEC_TARGET fptype device_TrigThresholdUpper (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x         = evt[indices[2 + indices[0]]]; 
  fptype thresh    = p[indices[1]];
  fptype trigConst = p[indices[2]];
  fptype linConst  = p[indices[3]];

  trigConst *= (thresh - x);
  return threshCalc(trigConst, linConst);
}

EXEC_TARGET fptype device_TrigThresholdLower (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x         = evt[indices[2 + indices[0]]]; 
  fptype thresh    = p[indices[1]];
  fptype trigConst = p[indices[2]];
  fptype linConst  = p[indices[3]];

  trigConst *= (x - thresh); 
  return threshCalc(trigConst, linConst);
}

EXEC_TARGET fptype device_VerySpecialEpisodeTrigThresholdUpper (fptype* evt, fptype* p, unsigned int* indices) {
  // Annoying special case for use with Mikhail's efficiency function across the Dalitz plot 

  fptype x         = evt[indices[2 + indices[0] + 0]]; 
  fptype y         = evt[indices[2 + indices[0] + 1]]; 

  fptype thresh    = p[indices[1]];
  fptype trigConst = p[indices[2]];
  fptype linConst  = p[indices[3]];
  fptype z         = p[indices[4]] - x - y; 

  trigConst *= (thresh - z); 
  return threshCalc(trigConst, linConst);
}

EXEC_TARGET fptype device_VerySpecialEpisodeTrigThresholdLower (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x         = evt[indices[2 + indices[0] + 0]]; 
  fptype y         = evt[indices[2 + indices[0] + 1]]; 

  fptype thresh    = p[indices[1]];
  fptype trigConst = p[indices[2]];
  fptype linConst  = p[indices[3]];
  fptype z         = p[indices[4]] - x - y; 

  trigConst *= (z - thresh); 
  fptype ret = threshCalc(trigConst, linConst);

  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //printf("TrigThreshold: (%f - %f = %f) -> %f %f\n", z, thresh, trigConst, linConst, ret);

  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_TrigThresholdUpper = device_TrigThresholdUpper; 
MEM_DEVICE device_function_ptr ptr_to_TrigThresholdLower = device_TrigThresholdLower; 
MEM_DEVICE device_function_ptr ptr_to_VerySpecialEpisodeTrigThresholdUpper = device_VerySpecialEpisodeTrigThresholdUpper; 
MEM_DEVICE device_function_ptr ptr_to_VerySpecialEpisodeTrigThresholdLower = device_VerySpecialEpisodeTrigThresholdLower; 


__host__ TrigThresholdPdf::TrigThresholdPdf (std::string n, Variable* _x, Variable* thresh, Variable* trigConst, Variable* linConst, bool upper) 
  : GooPdf(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(thresh));
  pindices.push_back(registerParameter(trigConst));
  pindices.push_back(registerParameter(linConst)); 

  if (upper) MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_TrigThresholdUpper, sizeof(void*), 0, cudaMemcpyDeviceToHost);
  else MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_TrigThresholdLower, sizeof(void*), 0, cudaMemcpyDeviceToHost);
  initialise(pindices); 
}

  __host__ TrigThresholdPdf::TrigThresholdPdf (std::string n, Variable* _x, Variable* _y, Variable* thresh, Variable* trigConst, Variable* linConst, Variable* massConstant, bool upper) 
  : GooPdf(0, n) 
{
  registerObservable(_x);
  registerObservable(_y);

  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(thresh));
  pindices.push_back(registerParameter(trigConst));
  pindices.push_back(registerParameter(linConst)); 
  pindices.push_back(registerParameter(massConstant)); 

  if (upper) MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_VerySpecialEpisodeTrigThresholdUpper, sizeof(void*), 0, cudaMemcpyDeviceToHost);
  else MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_VerySpecialEpisodeTrigThresholdLower, sizeof(void*), 0, cudaMemcpyDeviceToHost);
  initialise(pindices); 
}
