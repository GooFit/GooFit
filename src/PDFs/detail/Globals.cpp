#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/detail/Globals.h>

namespace GooFit {

void *host_fcn_ptr = nullptr;

std::map<void *, int> functionAddressToDeviceIndexMap;

fptype *dev_event_array;

fptype host_parameters[GOOFIT_MAXPAR];
fptype host_constants[GOOFIT_MAXPAR];
fptype host_observables[GOOFIT_MAXPAR];
fptype host_normalizations[GOOFIT_MAXPAR];

int host_callnumber     = 0;
int totalParameters     = 0;
int totalConstants      = 0;
int totalObservables    = 0;
int totalNormalizations = 0;

void *host_function_table[GOOFIT_MAXFUNC];
unsigned int num_device_functions = 0;

} // namespace GooFit
