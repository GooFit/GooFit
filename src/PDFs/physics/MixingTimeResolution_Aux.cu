#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

MixingTimeResolution::MixingTimeResolution()
    : GooPdf("mixing_resolution") {}
MixingTimeResolution::~MixingTimeResolution() = default;

void MixingTimeResolution::initIndex(void *dev_fcn_ptr) { resFunctionIdx = GooPdf::findFunctionIdx(dev_fcn_ptr); }

void MixingTimeResolution::createParameters(GooFit::PdfBase *) {}

} // namespace GooFit
