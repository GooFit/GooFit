#pragma once

#include <vector>

#include <mcbooster/GContainers.h>

#include <goofit/detail/Complex.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>

namespace GooFit {

class NormEvents_4Body_Weighted final : public NormEvents_4Body_Base {
  public:
    NormEvents_4Body_Weighted()                                                   = delete;
    NormEvents_4Body_Weighted(const NormEvents_4Body_Weighted &copyMe)            = default;
    NormEvents_4Body_Weighted(NormEvents_4Body_Weighted &&moveMe)                 = default;
    virtual ~NormEvents_4Body_Weighted() override                                 = default;
    NormEvents_4Body_Weighted &operator=(const NormEvents_4Body_Weighted &copyMe) = default;
    NormEvents_4Body_Weighted &operator=(NormEvents_4Body_Weighted &&moveMe)      = default;

    NormEvents_4Body_Weighted(mcbooster::RealVector_h m12,
                              mcbooster::RealVector_h m34,
                              mcbooster::RealVector_h cos12,
                              mcbooster::RealVector_h cos34,
                              mcbooster::RealVector_h phi,
                              mcbooster::RealVector_h dtime,
                              mcbooster::RealVector_h sigma,
                              mcbooster::RealVector_h weights);

    __host__ virtual fptype computeNorm_TD(bool noCachedNormValuesToCompute,
                                           const MixingTimeResolution *const resolution,
                                           fptype tau,
                                           fptype xmixing,
                                           fptype ymixing,
                                           unsigned int dalitzID,
                                           bool spinsCalculated,
                                           const std::vector<bool> &lineshapeChanged,
                                           const std::vector<unsigned int> &sfFunctionIndices,
                                           const std::vector<unsigned int> &lsFunctionIndices) override;

  private:
    // Store normalization events.
    mcbooster::RealVector_h _norm_M12_h;
    mcbooster::RealVector_h _norm_M34_h;
    mcbooster::RealVector_h _norm_CosTheta12_h;
    mcbooster::RealVector_h _norm_CosTheta34_h;
    mcbooster::RealVector_h _norm_phi_h;
    mcbooster::RealVector_h _norm_dtime_h;
    // TODO: Do we actually need this for the integral? Or should it drop out?
    mcbooster::RealVector_h _norm_sigma_h;
    // Store initial weight for normalization events.
    mcbooster::RealVector_h _norm_weight_h;
    // Store spin and lineshape values for normalization.
    mcbooster::RealVector_h _norm_SF_h;
    mcbooster::mc_host_vector<fpcomplex> _norm_LS_h;
};

} // namespace GooFit
