#include "MixingTimeResolution_Aux.hh"
#include "EngineCore.hh" 

MixingTimeResolution::MixingTimeResolution () {}
MixingTimeResolution::~MixingTimeResolution () {}

void MixingTimeResolution::initIndex (void* dev_fcn_ptr) {
  resFunctionIdx = EngineCore::findFunctionIdx(dev_fcn_ptr); 
}
