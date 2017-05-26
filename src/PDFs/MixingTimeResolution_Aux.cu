#include "goofit/PDFs/MixingTimeResolution_Aux.h"
#include "goofit/PDFs/GooPdf.h"

namespace GooFit {


MixingTimeResolution::MixingTimeResolution() {}
MixingTimeResolution::~MixingTimeResolution() {}

void MixingTimeResolution::initIndex(void* dev_fcn_ptr) {
    resFunctionIdx = GooPdf::findFunctionIdx(dev_fcn_ptr);
}
} // namespace GooFit

