#include <goofit/GlobalCudaDefines.h>

// TODO: These aren't the same as GooFit 2.0
// #define mPiPlus 0.13957018
#define mPiPlus 0.139570
#define mKPlus 0.493677
// #define mEta 0.54751
#define mEta 0.547862
#define mEtap 0.96778

namespace GooFit {

// Form factors as in pdg http://pdg.lbl.gov/2012/reviews/rpp2012-rev-dalitz-analysis-formalizm.pdf
// Does not account for L>2

// TODO: Non-macros should never be named with all caps!

__device__ fptype BL_PRIME(fptype z2, fptype z02, int L);
__device__ fptype BL(fptype z2, int L);
__device__ fptype BL2(fptype z2, int L);

} // namespace GooFit
