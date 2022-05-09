#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/detail/Globals.h>

namespace GooFit {

void *host_fcn_ptr = nullptr;

std::map<void *, int> functionAddressToDeviceIndexMap;

fptype *dev_event_array;

int host_callnumber = 0;

} // namespace GooFit
