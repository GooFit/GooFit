#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/detail/Globals.h>

namespace GooFit {

void *host_fcn_ptr = nullptr;

std::map<void *, int> functionAddressToDeviceIndexMap;

} // namespace GooFit
