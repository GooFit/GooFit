#pragma once

#include <vector>

#include <mcbooster/GContainers.h>

#include <goofit/detail/Complex.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>

namespace GooFit {

class NormEvents_4Body_HostCached final : public NormEvents_4Body_Base {
  public:
    NormEvents_4Body_HostCached()                                                     = delete;
    NormEvents_4Body_HostCached(const NormEvents_4Body_HostCached &copyMe)            = default;
    NormEvents_4Body_HostCached(NormEvents_4Body_HostCached &&moveMe)                 = default;
    virtual ~NormEvents_4Body_HostCached() override                                   = default;
    NormEvents_4Body_HostCached &operator=(const NormEvents_4Body_HostCached &copyMe) = default;
    NormEvents_4Body_HostCached &operator=(NormEvents_4Body_HostCached &&moveMe)      = default;

    static std::vector<NormEvents_4Body_Base *>
    buildBatches(const std::vector<long> &normSeeds,
                 unsigned int numNormEventsToGenPerBatch,
                 const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses);

    __host__ NormEvents_4Body_HostCached(long normSeed,
                                         unsigned int numNormEventsToGenPerBatch,
                                         const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses);

    __host__ virtual fptype computeNorm_TD(bool noCachedNormValuesToCompute,
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
    mcbooster::RealVector_h _norm_M12_h;
    mcbooster::RealVector_h _norm_M34_h;
    mcbooster::RealVector_h _norm_CosTheta12_h;
    mcbooster::RealVector_h _norm_CosTheta34_h;
    mcbooster::RealVector_h _norm_phi_h;
    // store spin and lineshape values for normalization
    mcbooster::RealVector_h _norm_SF_h;
    mcbooster::mc_host_vector<fpcomplex> _norm_LS_h;
};

} // end namespace GooFit
