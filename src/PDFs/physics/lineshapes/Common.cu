#include <goofit/GlobalCudaDefines.h>

#include "Common.h"

namespace GooFit {

// Form factors as in pdg http://pdg.lbl.gov/2012/reviews/rpp2012-rev-dalitz-analysis-formalizm.pdf
__device__ auto BL_PRIME(fptype z2, fptype z02, int L) -> fptype {
    if(0 == L)
        return 1.0;
    else if(1 == L)
        return (1 + z02) / (1 + z2);
    else if(2 == L)
        return (z02 * z02 + 3 * z02 + 9) / (z2 * z2 + 3 * z2 + 9);
    else {
        printf("ERROR! Oribtal > 2 not supported!\n");
        return 0;
    }

    // Spin 3 and up not accounted for.
}

__device__ auto BL(fptype z2, int L) -> fptype {
    if(0 == L)
        return 1.0;
    else if(1 == L)
        return 2 * z2 / (1 + z2);
    else if(2 == L)
        return (13 * z2 * z2) / (z2 * z2 + 3 * z2 + 9);
    else {
        printf("ERROR! Oribtal > 2 not supported!\n");
        return 0;
    }

    // Spin 3 and up not accounted for.
}

__device__ auto BL2(fptype z2, int L) -> fptype {
    if(0 == L)
        return 1.0;
    else if(1 == L)
        return 1.0 / (1 + z2);
    else if(2 == L)
        return 1.0 / (z2 * z2 + 3 * z2 + 9);
    else {
        printf("ERROR! Oribtal > 2 not supported!\n");
        return 0;
    }

    // Spin 3 and up not accounted for.
}

} // namespace GooFit
