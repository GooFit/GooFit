#include <Eigen/Core>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/detail/Complex.h>

#define NCHANNELS 5

namespace GooFit {
__device__ fpcomplex phsp_twoBody(fptype s, fptype m0, fptype m1);

__device__ fpcomplex phsp_fourPi(fptype s);

__device__ Eigen::Array<fpcomplex, NCHANNELS, NCHANNELS>
getPropagator(const Eigen::Array<fptype, NCHANNELS, NCHANNELS> &kMatrix,
              const Eigen::Matrix<fpcomplex, 5, 1> &phaseSpace,
              fptype adlerTerm);

} // namespace GooFit
