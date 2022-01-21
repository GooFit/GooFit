#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/SmartVector.h>

#include <map>

namespace GooFit {

struct ParameterContainer;

typedef fptype (*device_function_ptr)(fptype *, ParameterContainer &);

extern SmartVector<fptype> host_parameters;
extern SmartVector<fptype> host_constants;
extern SmartVector<fptype> host_observables;
extern SmartVector<fptype> host_normalizations;
extern SmartVector<void *> host_function_table;

// 211217  mds additions
extern std::vector<std::string> host_function_name;
extern std::vector<std::string> host_function_pdf;
extern std::vector<std::string> host_parameter_name;
extern std::vector<std::string> host_observable_name;

extern __device__ fptype *d_parameters;
extern __device__ fptype *d_constants;
extern __device__ fptype *d_observables;
extern __device__ fptype *d_normalizations;
extern __device__ void **d_function_table;

__host__ void cleanup();
__device__ auto dev_powi(int base, int exp) -> int; // Used in SmoothHistogramPdf and InterHistPdf

// a couple constants
extern __constant__ fptype c_motherMass;
extern __constant__ fptype c_daug1Mass;
extern __constant__ fptype c_daug2Mass;
extern __constant__ fptype c_daug3Mass;
extern __constant__ fptype c_meson_radius;
extern __constant__ fptype c_mother_meson_radius;

extern __constant__ unsigned int c_totalEvents;

extern std::map<void *, int> functionAddressToDeviceIndexMap;
extern std::map<void *, std::string> functionPtrToNameMap;

extern void *host_fcn_ptr;
extern fptype *dev_event_array;

} // namespace GooFit
