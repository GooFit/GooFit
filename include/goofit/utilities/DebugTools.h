#pragma once

#include <vector>

namespace GooFit {

  class DebugTools final {

  public:
    static __host__ std::vector<unsigned int> copyAmpIndicesToHost();

  };

} // end namespace GooFit
