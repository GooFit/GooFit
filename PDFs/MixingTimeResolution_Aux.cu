#include "MixingTimeResolution_Aux.hh"
#include "ThrustPdfFunctor.hh" 

MixingTimeResolution::MixingTimeResolution () {}
MixingTimeResolution::~MixingTimeResolution () {}

void MixingTimeResolution::initIndex (void* dev_fcn_ptr) {
  resFunctionIdx = ThrustPdfFunctor::findFunctionIdx(dev_fcn_ptr); 
}
