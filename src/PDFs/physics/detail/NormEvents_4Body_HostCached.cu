#include <numeric>
#include <vector>
#include <tuple>

#include <thrust/copy.h>

#include <mcbooster/GContainers.h>
#include <mcbooster/GTypes.h>

#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_HostCached.h>
#include <goofit/MathUtils.h>

namespace GooFit {

  __host__ NormEvents_4Body_HostCached::NormEvents_4Body_HostCached(
							   const std::vector<long>& normSeeds,
							   unsigned int numNormEventsToGenPerBatch,
							   const std::vector<mcbooster::GReal_t>& motherAndDaughterMasses)
  {
    _numBatches = normSeeds.size();

    _norm_M12_h = mcbooster::RealVector_h(0);
    _norm_M34_h = mcbooster::RealVector_h(0);
    _norm_CosTheta12_h = mcbooster::RealVector_h(0);
    _norm_CosTheta34_h = mcbooster::RealVector_h(0);
    _norm_phi_h = mcbooster::RealVector_h(0);

    _norm_SF_h = mcbooster::RealVector_h(0);
    _norm_LS_h = mcbooster::mc_host_vector<fpcomplex>(0);

    _totNumAccNormEvents = 0;
    _numAccPerBatch = std::vector<unsigned int>(_numBatches);
    
    // generate norm events on device in batches then save to host
    for (int b = 0; b < _numBatches; b++)
    {
      mcbooster::RealVector_d temp_norm_M12_d(0);
      mcbooster::RealVector_d temp_norm_M34_d(0);
      mcbooster::RealVector_d temp_norm_CosTheta12_d(0);
      mcbooster::RealVector_d temp_norm_CosTheta34_d(0);
      mcbooster::RealVector_d temp_norm_phi_d(0);

      unsigned int nAccThisBatch = NormEvents_4Body_Base::generate4BodyNormEvents(normSeeds[b],
										  numNormEventsToGenPerBatch,
										  motherAndDaughterMasses,
										  temp_norm_M12_d,
										  temp_norm_M34_d,
										  temp_norm_CosTheta12_d,
										  temp_norm_CosTheta34_d,
										  temp_norm_phi_d);
         
      _totNumAccNormEvents += nAccThisBatch;
      _numAccPerBatch[b] = nAccThisBatch;

      _norm_M12_h.resize(_totNumAccNormEvents);
      _norm_M34_h.resize(_totNumAccNormEvents);
      _norm_CosTheta12_h.resize(_totNumAccNormEvents);
      _norm_CosTheta34_h.resize(_totNumAccNormEvents);
      _norm_phi_h.resize(_totNumAccNormEvents);
     
      thrust::copy_n(temp_norm_M12_d.cbegin(),
		     nAccThisBatch, 
		     _norm_M12_h.begin()+getNumAccBeforeThisBatch(b));
      thrust::copy_n(temp_norm_M34_d.cbegin(), 
		     nAccThisBatch, 
		     _norm_M34_h.begin()+getNumAccBeforeThisBatch(b));
      thrust::copy_n(temp_norm_CosTheta12_d.cbegin(),
		     nAccThisBatch, 
		     _norm_CosTheta12_h.begin()+getNumAccBeforeThisBatch(b));
      thrust::copy_n(temp_norm_CosTheta34_d.cbegin(), 
		     nAccThisBatch, 
		     _norm_CosTheta34_h.begin()+getNumAccBeforeThisBatch(b));
      thrust::copy_n(temp_norm_phi_d.cbegin(), 
		     nAccThisBatch,
		     _norm_phi_h.begin()+getNumAccBeforeThisBatch(b));      

      GOOFIT_INFO("# of accepted MC events used for normalization [batch #{}]: {}", b, _numAccPerBatch[b]);
    }

    GOOFIT_INFO("Total # of accepted MC events used for normalization: {} ({} batches)", getNumAccNormEvents(), getNumBatches());
  }


  unsigned int NormEvents_4Body_HostCached::getNumAccBeforeThisBatch(unsigned int batchNum) const
  {
    return std::accumulate(
			   _numAccPerBatch.cbegin(), 
			   _numAccPerBatch.cbegin()+batchNum,
			   0);
  }


  __host__ fptype NormEvents_4Body_HostCached::computeNorm_TD(
							      bool noCachedNormValuesToCompute,
							      const MixingTimeResolution* const resolution,
							      fptype tau,
							      fptype xmixing,
							      fptype ymixing,
							      unsigned int dalitzId,					     
							      bool spinsCalculated,
							      const std::vector<bool>& lineshapeChanged,
							      const std::vector<unsigned int>& sfFunctionIndices,
							      const std::vector<unsigned int>& lsFunctionIndices)
  {
    // resize if needed
    if (!noCachedNormValuesToCompute)
    {
      unsigned int numSFCacheEntries = _totNumAccNormEvents * sfFunctionIndices.size();
      if (_norm_SF_h.size() != numSFCacheEntries)
      {
	_norm_SF_h.resize(numSFCacheEntries);
      }

      unsigned int numLSCacheEntries = _totNumAccNormEvents * lsFunctionIndices.size();
      if (_norm_LS_h.size() != numLSCacheEntries)
      {
	_norm_LS_h.resize(numLSCacheEntries);
      }
    }

    // do computations on device in batches
    std::vector<fptype> normResults(_numBatches);
    for (int b = 0; b < getNumBatches(); b++)
    {
      unsigned int numAccThisBatch = _numAccPerBatch[b];

      // copy previously cached sf values for batch to device
      unsigned int numSF = sfFunctionIndices.size();
      unsigned int numSFEntriesThisBatch = numAccThisBatch * numSF;
      unsigned int sfBeginOffset = getNumAccBeforeThisBatch(b) * numSF;
      mcbooster::RealVector_d sf_batchResult_d(numSFEntriesThisBatch);
      thrust::copy_n(_norm_SF_h.cbegin()+sfBeginOffset,
                     numSFEntriesThisBatch,
                     sf_batchResult_d.begin());

      // copy previously cached ls values for batch to device
      unsigned int numLS = lsFunctionIndices.size();
      unsigned int numLSEntriesThisBatch = numAccThisBatch * numLS;
      unsigned int lsBeginOffset = getNumAccBeforeThisBatch(b) * numLS;
      mcbooster::mc_device_vector<fpcomplex> ls_batchResult_d(numLSEntriesThisBatch);
      thrust::copy_n(_norm_LS_h.cbegin()+lsBeginOffset,
                     numLSEntriesThisBatch,
                     ls_batchResult_d.begin());

      if (!noCachedNormValuesToCompute)
      {
	// copy batch of norm events to device
	unsigned int beginOffset = getNumAccBeforeThisBatch(b);
	unsigned int endOffset = beginOffset + numAccThisBatch;
	mcbooster::RealVector_d batchNormM12_d(_norm_M12_h.cbegin()+beginOffset,
					       _norm_M12_h.cbegin()+endOffset);
	mcbooster::RealVector_d batchNormM34_d(_norm_M34_h.cbegin()+beginOffset,
					       _norm_M34_h.cbegin()+endOffset);
	mcbooster::RealVector_d batchNormCosTheta12_d(_norm_CosTheta12_h.cbegin()+beginOffset,
						      _norm_CosTheta12_h.cbegin()+endOffset);
	mcbooster::RealVector_d batchNormCosTheta34_d(_norm_CosTheta34_h.cbegin()+beginOffset,
						      _norm_CosTheta34_h.cbegin()+endOffset);
	mcbooster::RealVector_d batchNormPhi_d(_norm_phi_h.cbegin()+beginOffset,
					       _norm_phi_h.cbegin()+endOffset);

	// recompute cached sf values for batch (if needed)
	bool sfValuesUpdated = NormEvents_4Body_Base::computeCachedSFValuesForBatch_TD(
										       spinsCalculated,
										       dalitzId,
										       sfFunctionIndices,
										       numAccThisBatch,
										       batchNormM12_d,
										       batchNormM34_d,
										       batchNormCosTheta12_d,
										       batchNormCosTheta34_d,
										       batchNormPhi_d,
										       0,
										       sf_batchResult_d);
     
	// copy sf values for batch back to host if any of the values changed
	if (sfValuesUpdated)
	{
	  thrust::copy_n(sf_batchResult_d.cbegin(),
			 numSFEntriesThisBatch,
			 _norm_SF_h.begin()+sfBeginOffset);
	}

	// recompute cached ls values for batch (if needed)
	bool lsValuesUpdated = NormEvents_4Body_Base::computeCachedLSValuesForBatch_TD(
										       lineshapeChanged,
										       dalitzId,
										       lsFunctionIndices,
										       numAccThisBatch,
										       batchNormM12_d,
										       batchNormM34_d,
										       batchNormCosTheta12_d,
										       batchNormCosTheta34_d,
										       batchNormPhi_d,
										       0,
										       ls_batchResult_d);

	// copy ls values back to host if any of the values changed
	if (lsValuesUpdated)
	{
	  thrust::copy_n(ls_batchResult_d.cbegin(),
			 numLSEntriesThisBatch,
			 _norm_LS_h.begin()+lsBeginOffset);
	}
      } // end if computing cached values

      // do norm integral
      normResults[b] = NormEvents_4Body_Base::doNormIntegral_TD(
								resolution,
								tau,
								xmixing,
								ymixing,
								dalitzId,
								numAccThisBatch,
								sf_batchResult_d,
								ls_batchResult_d);
    } // end loop over batches
    
    fptype normResultsSum = MathUtils::doNeumaierSummation(normResults);
    
    return normResultsSum / _totNumAccNormEvents;
  }

} // end namespace GooFit
