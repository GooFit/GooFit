#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>

namespace GooFit {

MixingTimeResolution::MixingTimeResolution()  = default;
MixingTimeResolution::~MixingTimeResolution() = default;

void MixingTimeResolution::initIndex(void *dev_fcn_ptr) { resFunctionIdx = GooPdf::findFunctionIdx(dev_fcn_ptr); }
} // namespace GooFit
