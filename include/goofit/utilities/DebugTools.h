#pragma once

#include <vector>
#include <thrust/device_vector.h>

#include <goofit/detail/Complex.h>

namespace GooFit {

class DebugTools final {
  public:
    static std::vector<unsigned int> copyAmpIndicesToHost();

    static void printDeviceVecComplexVals(
        thrust::device_vector<fpcomplex>::const_iterator first,
        thrust::device_vector<fpcomplex>::const_iterator last,
        thrust::iterator_difference<thrust::device_vector<fpcomplex>::const_iterator>::type stride,
        const std::string &codeLoc,
        int printLimit,
        int printPrecision);
};

} // end namespace GooFit
