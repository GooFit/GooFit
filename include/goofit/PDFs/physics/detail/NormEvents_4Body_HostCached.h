#pragma once

#include <vector>

#include <mcbooster/GContainers.h>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>

namespace GooFit {

  class NormEvents_4Body_HostCached final : public NormEvents_4Body_Base {
  public:
    NormEvents_4Body_HostCached() = delete;
    NormEvents_4Body_HostCached(const NormEvents_4Body_HostCached& copyMe) = default;
    NormEvents_4Body_HostCached(NormEvents_4Body_HostCached&& moveMe) = default;
    ~NormEvents_4Body_HostCached() override = default;
    NormEvents_4Body_HostCached& operator=(const NormEvents_4Body_HostCached& copyMe) = default;
    NormEvents_4Body_HostCached& operator=(NormEvents_4Body_HostCached&& moveMe) = default;

    __host__ NormEvents_4Body_HostCached(
				const std::vector<long>& normSeeds,
				unsigned int numNormEventsToGenPerBatch,
				const std::vector<mcbooster::GReal_t>& motherAndDaughterMasses);

    __host__ virtual NormEvents_4Body_Batch getBatch(unsigned int batchNum) const override;

  protected:

  private:
    std::vector<unsigned int> _numAccPerBatch;
    // store normalization events
    mcbooster::RealVector_h _norm_M12_h;
    mcbooster::RealVector_h _norm_M34_h;
    mcbooster::RealVector_h _norm_CosTheta12_h;
    mcbooster::RealVector_h _norm_CosTheta34_h;
    mcbooster::RealVector_h _norm_phi_h;

    unsigned int getBatchOffset(unsigned int batchNum) const;
  };

} // end namespace GooFit
