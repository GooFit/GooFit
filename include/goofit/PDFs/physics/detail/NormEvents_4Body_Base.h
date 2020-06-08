#pragma once

#include <vector>

#include <mcbooster/GContainers.h>

#include <goofit/GlobalCudaDefines.h>

namespace GooFit {

  class NormEvents_4Body_Batch final {
  public:
    NormEvents_4Body_Batch(
			   const mcbooster::RealVector_d& norm_M12,
			   const mcbooster::RealVector_d& norm_M34,
			   const mcbooster::RealVector_d& norm_CosTheta12,
			   const mcbooster::RealVector_d& norm_CosTheta34,
			   const mcbooster::RealVector_d& norm_phi,
			   unsigned int numAccThisBatch,
			   unsigned int numAccBeforeThisBatch)
      : _NORM_M12(norm_M12),
      _NORM_M34(norm_M34),
      _NORM_COSTHETA12(norm_CosTheta12),
      _NORM_COSTHETA34(norm_CosTheta34),
      _NORM_PHI(norm_phi),
      _NUM_ACC_THIS_BATCH(numAccThisBatch),
      _NUM_ACC_BEFORE_THIS_BATCH(numAccBeforeThisBatch) {}

    NormEvents_4Body_Batch() = default;
    NormEvents_4Body_Batch(const NormEvents_4Body_Batch& copyMe) = default;
    NormEvents_4Body_Batch(NormEvents_4Body_Batch&& moveMe) = default;
    ~NormEvents_4Body_Batch() = default;
    NormEvents_4Body_Batch& operator=(const NormEvents_4Body_Batch& copyMe) = default;
    NormEvents_4Body_Batch& operator=(NormEvents_4Body_Batch&& moveMe) = default;

    const mcbooster::RealVector_d _NORM_M12;
    const mcbooster::RealVector_d _NORM_M34;
    const mcbooster::RealVector_d _NORM_COSTHETA12;
    const mcbooster::RealVector_d _NORM_COSTHETA34;
    const mcbooster::RealVector_d _NORM_PHI;
    const unsigned int _NUM_ACC_THIS_BATCH;
    const unsigned int _NUM_ACC_BEFORE_THIS_BATCH;
  protected:

  private:
  };


  class NormEvents_4Body_Base {
  public:
    NormEvents_4Body_Base() = default;
    NormEvents_4Body_Base(const NormEvents_4Body_Base& copyMe) = default;
    NormEvents_4Body_Base(NormEvents_4Body_Base&& moveMe) = default;
    virtual ~NormEvents_4Body_Base() = default;
    NormEvents_4Body_Base& operator=(const NormEvents_4Body_Base& copyMe) = default; 
    NormEvents_4Body_Base& operator=(NormEvents_4Body_Base&& moveMe) = default;

    int getNumAccNormEvents() const { return _totNumAccNormEvents; } 

    int getNumBatches() const { return _numBatches; }

    __host__ virtual NormEvents_4Body_Batch getBatch(unsigned int batchNum) const = 0;

  protected:
    unsigned int _totNumAccNormEvents;
    unsigned int _numBatches;

    __host__ static unsigned int generate4BodyNormEvents(
							 long normSeed,
							 unsigned int numNormEventsToGen,
							 const std::vector<mcbooster::GReal_t>& motherAndDaughterMasses,
							 mcbooster::RealVector_d& norm_M12,
							 mcbooster::RealVector_d& norm_M34,
							 mcbooster::RealVector_d& norm_CosTheta12,
							 mcbooster::RealVector_d& norm_CosTheta34,
							 mcbooster::RealVector_d& norm_phi);

  private:
  };

} // end namespace GooFit
