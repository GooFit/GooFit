#pragma once

#include <vector>

#include <mcbooster/GContainers.h>

#include <goofit/detail/Complex.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>

namespace GooFit {

class NormEvents_4Body_WeightedDevice final : public NormEvents_4Body_Base {
  public:
    NormEvents_4Body_WeightedDevice()                                                         = delete;
    NormEvents_4Body_WeightedDevice(const NormEvents_4Body_WeightedDevice &copyMe)            = default;
    NormEvents_4Body_WeightedDevice(NormEvents_4Body_WeightedDevice &&moveMe)                 = default;
    virtual ~NormEvents_4Body_WeightedDevice() override                                       = default;
    NormEvents_4Body_WeightedDevice &operator=(const NormEvents_4Body_WeightedDevice &copyMe) = default;
    NormEvents_4Body_WeightedDevice &operator=(NormEvents_4Body_WeightedDevice &&moveMe)      = default;

    NormEvents_4Body_WeightedDevice(mcbooster::RealVector_d m12,
                                    mcbooster::RealVector_d m34,
                                    mcbooster::RealVector_d cos12,
                                    mcbooster::RealVector_d cos34,
                                    mcbooster::RealVector_d phi,
                                    mcbooster::RealVector_d dtime,
                                    mcbooster::RealVector_d sigma,
                                    mcbooster::RealVector_d weights);

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
    mcbooster::RealVector_d _norm_M12_d;
    mcbooster::RealVector_d _norm_M34_d;
    mcbooster::RealVector_d _norm_CosTheta12_d;
    mcbooster::RealVector_d _norm_CosTheta34_d;
    mcbooster::RealVector_d _norm_phi_d;
    mcbooster::RealVector_d _norm_dtime_d;
    // TODO: Do we actually need this for the integral? Or should it drop out?
    mcbooster::RealVector_d _norm_sigma_d;
    // Store initial weight for normalization events.
    mcbooster::RealVector_d _norm_weight_d;
    // Store spin and lineshape values for normalization.
    mcbooster::RealVector_d _norm_SF_d;
    mcbooster::mc_device_vector<fpcomplex> _norm_LS_d;
};

} // namespace GooFit
