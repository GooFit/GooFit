#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

#include <thrust/device_vector.h>

#include "goofit/utilities/DebugTools.h"
#include "goofit/GlobalCudaDefines.h"
#include "goofit/PDFs/physics/Amp4BodyGlobals.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

void DebugTools::printDeviceVecComplexVals(
    thrust::device_vector<fpcomplex>::const_iterator first,
    thrust::device_vector<fpcomplex>::const_iterator last,
    thrust::iterator_difference<thrust::device_vector<fpcomplex>::const_iterator>::type stride,
    const std::string &codeLoc,
    int printLimit,
    int printPrecision) {
    int numPrinted = 0;

    auto pStridedRange = strided_range<thrust::device_vector<fpcomplex>::const_iterator>(first, last, stride);

    for(auto pIter = pStridedRange.begin(); pIter != pStridedRange.end(); pIter = std::next(pIter)) {
        fpcomplex entry = *pIter;
        std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(printPrecision) << codeLoc << ": "
                  << entry.real() << ", " << entry.imag() << std::endl;

        numPrinted++;
        if(numPrinted > printLimit) {
            break;
        }
    }
    std::cout << std::endl;
}

std::vector<unsigned int> DebugTools::copyAmpIndicesToHost() {
    std::vector<unsigned int> hostAmpIndices(500);

    MEMCPY_FROM_SYMBOL(&(hostAmpIndices[0]), AmpIndices, 500 * sizeof(unsigned int), 0, cudaMemcpyDeviceToHost);

    return hostAmpIndices;
}

} // namespace GooFit
