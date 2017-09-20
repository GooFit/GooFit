#include "goofit/PDFs/physics/MixingTimeResolution_Aux.h"
#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

MixingTimeResolution::MixingTimeResolution() : GooPdf(NULL, "mixing_resolution") {}
MixingTimeResolution::~MixingTimeResolution() = default;

void MixingTimeResolution::initIndex(void *dev_fcn_ptr) { resFunctionIdx = GooPdf::findFunctionIdx(dev_fcn_ptr); }

void MixingTimeResolution::recursiveSetIndices() {
}
} // namespace GooFit
