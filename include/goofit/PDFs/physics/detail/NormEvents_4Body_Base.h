#pragma once

#include <vector>

#include <mcbooster/GContainers.h>

#include <goofit/detail/Complex.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

class NormEvents_4Body_Base {
  public:
    NormEvents_4Body_Base()                                               = default;
    NormEvents_4Body_Base(const NormEvents_4Body_Base &copyMe)            = default;
    NormEvents_4Body_Base(NormEvents_4Body_Base &&moveMe)                 = default;
    virtual ~NormEvents_4Body_Base()                                      = default;
    NormEvents_4Body_Base &operator=(const NormEvents_4Body_Base &copyMe) = default;
    NormEvents_4Body_Base &operator=(NormEvents_4Body_Base &&moveMe)      = default;

    int getNumAccNormEvents() const { return _totNumAccNormEvents; }

    __host__ virtual fptype computeNorm_TD(bool noCachedNormValuesToCompute,
                                           const MixingTimeResolution *const resolution,
                                           fptype tau,
                                           fptype xmixing,
                                           fptype ymixing,
                                           unsigned int dalitzId,
                                           bool spinsCalculated,
                                           const std::vector<bool> &lineshapeChanged,
                                           const std::vector<unsigned int> &sfFunctionIndices,
                                           const std::vector<unsigned int> &lsFunctionIndices)
        = 0;

  protected:
    __host__ static unsigned int generate4BodyNormEvents(long normSeed,
                                                         unsigned int numNormEventsToGen,
                                                         const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses,
                                                         mcbooster::RealVector_d &norm_M12,
                                                         mcbooster::RealVector_d &norm_M34,
                                                         mcbooster::RealVector_d &norm_CosTheta12,
                                                         mcbooster::RealVector_d &norm_CosTheta34,
                                                         mcbooster::RealVector_d &norm_phi);

    __host__ static bool computeCachedLSValuesForBatch_TD(const std::vector<bool> &lineshapeChanged,
                                                          unsigned int dalitzId,
                                                          const std::vector<unsigned int> &lsFunctionIndices,
                                                          unsigned int numAccThisBatch,
                                                          const mcbooster::RealVector_d &batchNormM12_d,
                                                          const mcbooster::RealVector_d &batchNormM34_d,
                                                          const mcbooster::RealVector_d &batchNormCosTheta12_d,
                                                          const mcbooster::RealVector_d &batchNormCosTheta34_d,
                                                          const mcbooster::RealVector_d &batchNormPhi_d,
                                                          unsigned int resultOffset,
                                                          mcbooster::mc_device_vector<fpcomplex> &ls_batchResult_d);

    __host__ static bool computeCachedSFValuesForBatch_TD(bool spinsCalculated,
                                                          unsigned int dalitzId,
                                                          const std::vector<unsigned int> &sfFunctionIndices,
                                                          unsigned int numAccThisBatch,
                                                          const mcbooster::RealVector_d &batchNormM12_d,
                                                          const mcbooster::RealVector_d &batchNormM34_d,
                                                          const mcbooster::RealVector_d &batchNormCosTheta12_d,
                                                          const mcbooster::RealVector_d &batchNormCosTheta34_d,
                                                          const mcbooster::RealVector_d &batchNormPhi_d,
                                                          unsigned int resultOffset,
                                                          mcbooster::RealVector_d &sf_batchResult_d);

    __host__ static fptype doNormIntegral_TD(const MixingTimeResolution *const resolution,
                                             fptype tau,
                                             fptype xmixing,
                                             fptype ymixing,
                                             unsigned int dalitzId,
                                             unsigned int numAccThisBatch,
                                             mcbooster::RealVector_d &batchSF_d,
                                             mcbooster::mc_device_vector<fpcomplex> &batchLS_d);

    unsigned int _totNumAccNormEvents;

  private:
};

} // end namespace GooFit
