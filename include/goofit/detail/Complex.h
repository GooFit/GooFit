#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <thrust/complex.h>

namespace GooFit {

using fpcomplex = thrust::complex<fptype>;
    
namespace literals {
namespace complex_literals {

inline fpcomplex operator"" _i (unsigned long long d) {
    return fpcomplex(0.0, d);
}

inline fpcomplex operator"" _i (long double d) {
    return fpcomplex(0.0, d);
}

}
}

} // GooFit
