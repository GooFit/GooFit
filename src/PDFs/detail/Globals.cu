#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/detail/Globals.h>

namespace GooFit {

// Device-side, translation-unit constrained.  These were constant, removing const.
// The reason is that this will make it much more difficult to fetch memory, since
// it has less memory to work with limiting the constant amount.
SmartVector<fptype> host_parameters{"d_parameters"};
SmartVector<fptype> host_constants{"d_constants"};
SmartVector<fptype> host_observables{"d_observables"};
SmartVector<fptype> host_normalizations{"d_normalizations"};
SmartVector<void *> host_function_table{"d_function_table"};

///  mds 211217
std::vector<std::string> host_function_name;
std::vector<std::string> host_function_pdf;
std::vector<std::string> host_parameter_name;
std::vector<std::string> host_observable_name;
///

__device__ fptype *d_parameters;
__device__ fptype *d_constants;
__device__ fptype *d_observables;
__device__ fptype *d_normalizations;
__device__ void **d_function_table;

__constant__ unsigned int c_totalEvents;
__constant__ fptype c_motherMass;
__constant__ fptype c_daug1Mass;
__constant__ fptype c_daug2Mass;
__constant__ fptype c_daug3Mass;
__constant__ fptype c_meson_radius;
__constant__ fptype c_mother_meson_radius;

/// Clear all device memory (call before exit!)
__host__ void cleanup() {
    host_parameters.clear_device();
    host_constants.clear_device();
    host_observables.clear_device();
    host_normalizations.clear_device();
    host_function_table.clear_device();
}

__device__ auto dev_powi(int base, int exp) -> int {
    int ret = 1;

    for(int i = 0; i < exp; ++i)
        ret *= base;

    return ret;
}

} // namespace GooFit
