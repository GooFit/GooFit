#include <numeric>
#include <vector>
#include <tuple>

#include <thrust/copy.h>

#include <mcbooster/GContainers.h>
#include <mcbooster/GTypes.h>

#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_HostCached.h>

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

    _totNumAccNormEvents = 0;
    _numAccPerBatch = std::vector<unsigned int>(_numBatches);
    
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
		     _norm_M12_h.begin()+getBatchOffset(b));
      thrust::copy_n(temp_norm_M34_d.cbegin(), 
		     nAccThisBatch, 
		     _norm_M34_h.begin()+getBatchOffset(b));
      thrust::copy_n(temp_norm_CosTheta12_d.cbegin(),
		     nAccThisBatch, 
		     _norm_CosTheta12_h.begin()+getBatchOffset(b));
      thrust::copy_n(temp_norm_CosTheta34_d.cbegin(), 
		     nAccThisBatch, 
		     _norm_CosTheta34_h.begin()+getBatchOffset(b));
      thrust::copy_n(temp_norm_phi_d.cbegin(), 
		     nAccThisBatch,
		     _norm_phi_h.begin()+getBatchOffset(b));      

      GOOFIT_INFO("# of accepted MC events used for normalization [batch #{}]: {}", b, _numAccPerBatch[b]);
    }

    GOOFIT_INFO("Total # of accepted MC events used for normalization: {} ({} batches)", getNumAccNormEvents(), getNumBatches());
  }


  unsigned int NormEvents_4Body_HostCached::getBatchOffset(unsigned int batchNum) const
  {
    return std::accumulate(
			   _numAccPerBatch.cbegin(), 
			   _numAccPerBatch.cbegin()+batchNum,
			   0);
  }


  __host__ NormEvents_4Body_Batch NormEvents_4Body_HostCached::getBatch(unsigned int batchNum) const
  { 
    if (batchNum >= _numBatches)
    {
      throw GooFit::GeneralError("Invalid batch number for NormEvents_4Body_HostCached object.");
    }

    unsigned int numAccThisBatch = _numAccPerBatch[batchNum];
    
    unsigned int beginOffset = getBatchOffset(batchNum);
    unsigned int endOffset = beginOffset + numAccThisBatch;

    mcbooster::RealVector_d norm_M12_d(_norm_M12_h.cbegin()+beginOffset,
				       _norm_M12_h.cbegin()+endOffset);
    mcbooster::RealVector_d norm_M34_d(_norm_M34_h.cbegin()+beginOffset,
				       _norm_M34_h.cbegin()+endOffset);
    mcbooster::RealVector_d norm_CosTheta12_d(_norm_CosTheta12_h.cbegin()+beginOffset,
					      _norm_CosTheta12_h.cbegin()+endOffset);
    mcbooster::RealVector_d norm_CosTheta34_d(_norm_CosTheta34_h.cbegin()+beginOffset,
					      _norm_CosTheta34_h.cbegin()+endOffset);
    mcbooster::RealVector_d norm_phi_d(_norm_phi_h.cbegin()+beginOffset,
				       _norm_phi_h.cbegin()+endOffset);

      
    return NormEvents_4Body_Batch(norm_M12_d,
				  norm_M34_d,
				  norm_CosTheta12_d,
				  norm_CosTheta34_d,
				  norm_phi_d,
				  numAccThisBatch,
				  beginOffset);
  }

} // end namespace GooFit
