#pragma once

#include <vector>

#include "goofit/GlobalCudaDefines.h"

namespace GooFit {

class MathUtils final {
  public:
    MathUtils()                                   = default;
    MathUtils(const MathUtils &copyMe)            = default;
    MathUtils(MathUtils &&moveMe)                 = default;
    virtual ~MathUtils()                          = default;
    MathUtils &operator=(const MathUtils &copyMe) = default;
    MathUtils &operator=(MathUtils &&moveMe)      = default;

    static fptype doNeumaierSummation(const std::vector<fptype> &vals);
};

} // end namespace GooFit
