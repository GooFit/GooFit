#pragma once

#include <mcbooster/GContainers.h>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>

namespace GooFit {

  class NormEvents_4Body_DeviceCached final : public NormEvents_4Body_Base {
  public:
    NormEvents_4Body_DeviceCached() = delete;
    NormEvents_4Body_DeviceCached(const NormEvents_4Body_DeviceCached& copyMe) = default;
    NormEvents_4Body_DeviceCached(NormEvents_4Body_DeviceCached&& moveMe) = default;
    ~NormEvents_4Body_DeviceCached() override = default;
    NormEvents_4Body_DeviceCached& operator=(const NormEvents_4Body_DeviceCached& copyMe) = default;
    NormEvents_4Body_DeviceCached& operator=(NormEvents_4Body_DeviceCached&& moveMe) = default;

    __host__ NormEvents_4Body_DeviceCached(
				  const std::vector<mcbooster::GReal_t>& motherAndDaughterMasses,
				  long normSeed,
				  unsigned int numNormEventsToGen);    

    __host__ virtual NormEvents_4Body_Batch getBatch(unsigned int batchNum) const override;

  protected:

  private:
    // store normalization events
    mcbooster::RealVector_d _norm_M12_d;
    mcbooster::RealVector_d _norm_M34_d;
    mcbooster::RealVector_d _norm_CosTheta12_d;
    mcbooster::RealVector_d _norm_CosTheta34_d;
    mcbooster::RealVector_d _norm_phi_d;
  };

} // end namespace GooFit
