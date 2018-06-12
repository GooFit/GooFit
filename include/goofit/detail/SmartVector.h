#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <vector>

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include "SmartVectorGPU.h"
#else
#include "SmartVectorCPU.h"
#endif
