#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/Version.h>

#include <map>

namespace GooFit {

struct ParameterContainer;

typedef fptype (*device_function_ptr)(fptype *, ParameterContainer &);

extern fptype host_parameters[]; // GOOFIT_MAXPAR
extern fptype host_constants[];
extern fptype host_observables[];
extern fptype host_normalizations[];

extern int totalParameters;
extern int totalConstants;
extern int totalObservables;
extern int totalNormalizations;

extern __device__ fptype d_parameters[GOOFIT_MAXPAR];
extern __device__ fptype d_constants[GOOFIT_MAXPAR];
extern __device__ fptype d_observables[GOOFIT_MAXPAR];
extern __device__ fptype d_normalizations[GOOFIT_MAXPAR];

extern __device__ void *device_function_table[GOOFIT_MAXFUNC];
extern void *host_function_table[GOOFIT_MAXFUNC];
extern unsigned int num_device_functions;

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
