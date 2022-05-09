#pragma once

#include <mcbooster/GContainers.h>

#include <goofit/detail/Complex.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

class NormEvents_4Body_DeviceCached final : public NormEvents_4Body_Base {
  public:
    NormEvents_4Body_DeviceCached()                                                       = delete;
    NormEvents_4Body_DeviceCached(const NormEvents_4Body_DeviceCached &copyMe)            = default;
    NormEvents_4Body_DeviceCached(NormEvents_4Body_DeviceCached &&moveMe)                 = default;
    virtual ~NormEvents_4Body_DeviceCached() override                                     = default;
    NormEvents_4Body_DeviceCached &operator=(const NormEvents_4Body_DeviceCached &copyMe) = default;
    NormEvents_4Body_DeviceCached &operator=(NormEvents_4Body_DeviceCached &&moveMe)      = default;

    static std::vector<NormEvents_4Body_Base *>
    buildBatches(const std::vector<long> &normSeeds,
                 unsigned int numNormEventsToGenPerBatch,
                 const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses);

    __host__ NormEvents_4Body_DeviceCached(const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses,
                                           long normSeed,
                                           unsigned int numNormEventsToGen);

    __host__ fptype virtual computeNorm_TD(bool noCachedNormValuesToCompute,
                                           const MixingTimeResolution *const resolution,
                                           fptype tau,
                                           fptype xmixing,
                                           fptype ymixing,
                                           unsigned int dalitzId,
                                           bool spinsCalculated,
                                           const std::vector<bool> &lineshapeChanged,
                                           const std::vector<unsigned int> &sfFunctionIndices,
                                           const std::vector<unsigned int> &lsFunctionIndices) override;

  protected:
  private:
    // store normalization events
    mcbooster::RealVector_d _norm_M12_d;
    mcbooster::RealVector_d _norm_M34_d;
    mcbooster::RealVector_d _norm_CosTheta12_d;
    mcbooster::RealVector_d _norm_CosTheta34_d;
    mcbooster::RealVector_d _norm_phi_d;
    // store spin and lineshape values for normalization
    mcbooster::RealVector_d _norm_SF_d;
    mcbooster::mc_device_vector<fpcomplex> _norm_LS_d;
};

} // end namespace GooFit
