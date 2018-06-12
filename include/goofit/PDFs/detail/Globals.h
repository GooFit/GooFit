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

extern __device__ fptype *d_parameters;
extern __device__ fptype *d_constants;
extern __device__ fptype *d_observables;
extern __device__ fptype *d_normalizations;
extern __device__ void **d_function_table;

__host__ void cleanup();
__device__ int dev_powi(int base, int exp); // Used in SmoothHistogramPdf and InterHistPdf

// a couple constants
extern __constant__ fptype c_motherMass;
extern __constant__ fptype c_daug1Mass;
extern __constant__ fptype c_daug2Mass;
extern __constant__ fptype c_daug3Mass;
extern __constant__ fptype c_meson_radius;

extern __constant__ unsigned int c_totalEvents;

extern std::map<void *, int> functionAddressToDeviceIndexMap;

extern void *host_fcn_ptr;
extern fptype *dev_event_array;

} // namespace GooFit
