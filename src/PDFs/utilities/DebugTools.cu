#include <vector>

#include "goofit/utilities/DebugTools.h"
#include "goofit/GlobalCudaDefines.h"
#include "goofit/PDFs/physics/Amp4BodyGlobals.h"

namespace GooFit {

__host__ std::vector<unsigned int> DebugTools::copyAmpIndicesToHost()
{
  std::vector<unsigned int> hostAmpIndices(500);

  MEMCPY_FROM_SYMBOL(&(hostAmpIndices[0]),
		     AmpIndices,
		     500*sizeof(unsigned int),
		     0,
		     cudaMemcpyDeviceToHost);

  return hostAmpIndices;
}

} // end namespace Goofit
