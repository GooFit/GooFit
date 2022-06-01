#include <cstdlib>

#include "goofit/GlobalCudaDefines.h"
#include "goofit/MathUtils.h"

namespace GooFit {

/**
 * See https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
fptype MathUtils::doNeumaierSummation(const std::vector<fptype> &vals) {
    if(vals.size() == 1) {
        return vals[0];
    }

    fptype sum = 0.0;
    fptype c   = 0.0; // A running compensation for lost low-order bits.

    for(int i = 0; i < vals.size(); i++) {
        fptype t = sum + vals[i];
        if(std::abs(sum) >= std::abs(vals[i])) {
            c += (sum - t) + vals[i]; // If sum is bigger, low-order digits of input[i] are lost.
        } else {
            c += (vals[i] - t) + sum; // Else low-order digits of sum are lost.
        }
        sum = t;
    }

    return sum + c; // Correction only applied once in the very end.
}

} // end namespace GooFit
