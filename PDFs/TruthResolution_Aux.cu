#include "TruthResolution_Aux.hh" 

EXEC_TARGET fptype device_truth_resolution (fptype coshterm, fptype costerm, fptype sinhterm, fptype sinterm, 
					   fptype tau, fptype dtime, fptype xmixing, fptype ymixing, fptype /*sigma*/, 
					   fptype* /*p*/, unsigned int* /*indices*/) { 
  fptype ret = 0;
  dtime /= tau; 
  ret += coshterm*COSH(ymixing * dtime);
  ret += costerm*COS (xmixing * dtime);
  ret -= 2*sinhterm * SINH(ymixing * dtime);
  ret -= 2*sinterm * SIN (xmixing * dtime); // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B. 
  ret *= EXP(-dtime); 

  //cuPrintf("device_truth_resolution %f %f %f %f %f\n", coshterm, costerm, sinhterm, sinterm, dtime); 
  return ret; 
}

MEM_DEVICE device_resfunction_ptr ptr_to_truth = device_truth_resolution; 

TruthResolution::TruthResolution () 
  : MixingTimeResolution()
{
  MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_truth, sizeof(void*), 0, cudaMemcpyDeviceToHost);
  initIndex(); 
}
TruthResolution::~TruthResolution () {} 

fptype TruthResolution::normalisation (fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const {
  fptype timeIntegralOne = tau / (1 - ymixing*ymixing); 
  fptype timeIntegralTwo = tau / (1 + xmixing*xmixing);
  fptype timeIntegralThr = ymixing * timeIntegralOne;
  fptype timeIntegralFou = xmixing * timeIntegralTwo;
       
  fptype ret = timeIntegralOne * (di1 + di2);
  ret       += timeIntegralTwo * (di1 - di2);
  ret       -= 2*timeIntegralThr * di3;
  ret       -= 2*timeIntegralFou * di4;

  return ret; 
}
