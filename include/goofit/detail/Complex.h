#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <thrust/complex.h>

namespace GooFit {

using fpcomplex = thrust::complex<fptype>;

namespace literals {
namespace complex_literals {

inline auto operator"" _i(unsigned long long d) -> fpcomplex { return {0.0, static_cast<fptype>(d)}; }

inline auto operator"" _i(long double d) -> fpcomplex { return {0.0, static_cast<fptype>(d)}; }

} // namespace complex_literals
} // namespace literals

} // namespace GooFit
