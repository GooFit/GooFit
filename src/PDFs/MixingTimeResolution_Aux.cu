#include "MixingTimeResolution_Aux.hh"
#include "GooPdf.hh" 

MixingTimeResolution::MixingTimeResolution () {}
MixingTimeResolution::~MixingTimeResolution () {}

void MixingTimeResolution::initIndex (void* dev_fcn_ptr) {
  resFunctionIdx = GooPdf::findFunctionIdx(dev_fcn_ptr); 
}
