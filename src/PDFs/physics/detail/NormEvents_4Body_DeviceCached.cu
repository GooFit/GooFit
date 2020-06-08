#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_DeviceCached.h>

namespace GooFit {

  __host__ NormEvents_4Body_DeviceCached::NormEvents_4Body_DeviceCached(
									const std::vector<mcbooster::GReal_t>& motherAndDaughterMasses,
									long normSeed,
									unsigned int numNormEventsToGen)
  {
    _numBatches = 1;

    _norm_M12_d = mcbooster::RealVector_d(0);
    _norm_M34_d = mcbooster::RealVector_d(0);
    _norm_CosTheta12_d = mcbooster::RealVector_d(0);
    _norm_CosTheta34_d = mcbooster::RealVector_d(0);
    _norm_phi_d = mcbooster::RealVector_d(0);

    _norm_SF_d = mcbooster::RealVector_d(0);
    _norm_LS_d = mcbooster::mc_device_vector<fpcomplex>(0);

    _totNumAccNormEvents = NormEvents_4Body_Base::generate4BodyNormEvents(normSeed,
  									  numNormEventsToGen,
  									  motherAndDaughterMasses,
  									  _norm_M12_d,
  									  _norm_M34_d,
  									  _norm_CosTheta12_d,
  									  _norm_CosTheta34_d,
  									  _norm_phi_d);

    GOOFIT_INFO("# of accepted MC events used for normalization: {}", getNumAccNormEvents());
  }


  /* for testing, to compare results with normevents_4body_hostcached */
  /* __host__ NormEvents_4Body_DeviceCached::NormEvents_4Body_DeviceCached( */
  /* 									const std::vector<mcbooster::GReal_t>& motherAndDaughterMasses, */
  /* 									long normSeed, */
  /* 									unsigned int numNormEventsToGen) */
  /* { */
  /*   _numBatches = 1; */
  /*   _norm_M12_d = mcbooster::RealVector_d(0); */
  /*   _norm_M34_d = mcbooster::RealVector_d(0); */
  /*   _norm_CosTheta12_d = mcbooster::RealVector_d(0); */
  /*   _norm_CosTheta34_d = mcbooster::RealVector_d(0); */
  /*   _norm_phi_d = mcbooster::RealVector_d(0); */
  /*   _norm_SF_d = mcbooster::RealVector_d(0); */
  /*   _norm_LS_d = mcbooster::mc_device_vector<fpcomplex>(0); */

  /*   unsigned int numAccSoFar = 0; */
  /*   for (int b = 0; b < 3; b++) */
  /*   { */
  /*     mcbooster::RealVector_d temp_norm_M12_d(0); */
  /*     mcbooster::RealVector_d temp_norm_M34_d(0); */
  /*     mcbooster::RealVector_d temp_norm_CosTheta12_d(0); */
  /*     mcbooster::RealVector_d temp_norm_CosTheta34_d(0); */
  /*     mcbooster::RealVector_d temp_norm_phi_d(0); */

  /*     unsigned int nAccThisBatch = NormEvents_4Body_Base::generate4BodyNormEvents(normSeed+b, */
  /* 										  numNormEventsToGen, */
  /* 										  motherAndDaughterMasses, */
  /* 										  temp_norm_M12_d, */
  /* 										  temp_norm_M34_d, */
  /* 										  temp_norm_CosTheta12_d, */
  /* 										  temp_norm_CosTheta34_d, */
  /* 										  temp_norm_phi_d); */
  /*     unsigned int numAccBeforeThisBatch = numAccSoFar; */
  /*     numAccSoFar += nAccThisBatch; */

  /*     _norm_M12_d.resize(numAccSoFar); */
  /*     _norm_M34_d.resize(numAccSoFar); */
  /*     _norm_CosTheta12_d.resize(numAccSoFar); */
  /*     _norm_CosTheta34_d.resize(numAccSoFar); */
  /*     _norm_phi_d.resize(numAccSoFar); */

  /*     thrust::copy_n(temp_norm_M12_d.cbegin(), */
  /* 		     nAccThisBatch, */
  /* 		     _norm_M12_d.begin()+numAccBeforeThisBatch); */
  /*     thrust::copy_n(temp_norm_M34_d.cbegin(), */
  /*                    nAccThisBatch, */
  /*                    _norm_M34_d.begin()+numAccBeforeThisBatch); */
  /*     thrust::copy_n(temp_norm_CosTheta12_d.cbegin(), */
  /*                    nAccThisBatch, */
  /*                    _norm_CosTheta12_d.begin()+numAccBeforeThisBatch); */
  /*     thrust::copy_n(temp_norm_CosTheta34_d.cbegin(), */
  /*                    nAccThisBatch, */
  /*                    _norm_CosTheta34_d.begin()+numAccBeforeThisBatch); */
  /*     thrust::copy_n(temp_norm_phi_d.cbegin(), */
  /*                    nAccThisBatch, */
  /*                    _norm_phi_d.begin()+numAccBeforeThisBatch); */
  /*   } */

  /*   _totNumAccNormEvents = numAccSoFar; */
  /*   GOOFIT_INFO("# of accepted MC events used for normalization: {}", getNumAccNormEvents()); */
  /* } */


  __host__ fptype NormEvents_4Body_DeviceCached::computeNorm_TD(
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
    if (!noCachedNormValuesToCompute)
    {
      unsigned int numSFCacheEntries = _totNumAccNormEvents * sfFunctionIndices.size();
      if (_norm_SF_d.size() != numSFCacheEntries)
      {
	_norm_SF_d.resize(numSFCacheEntries);
      }

      // compute cached sf values (if needed)
      NormEvents_4Body_Base::computeCachedSFValuesForBatch_TD(
							      spinsCalculated,
							      dalitzId,
							      sfFunctionIndices,
							      _totNumAccNormEvents,
							      _norm_M12_d,
							      _norm_M34_d,
							      _norm_CosTheta12_d,
							      _norm_CosTheta34_d,
							      _norm_phi_d,
							      0,
							      _norm_SF_d);
      
      unsigned int numLSCacheEntries = _totNumAccNormEvents * lsFunctionIndices.size();
      if (_norm_LS_d.size() != numLSCacheEntries)
      {
	_norm_LS_d.resize(numLSCacheEntries);
      }

      // recompute cached ls values for batch (if needed)
      NormEvents_4Body_Base::computeCachedLSValuesForBatch_TD(
							      lineshapeChanged,
							      dalitzId,
							      lsFunctionIndices,
							      _totNumAccNormEvents,
							      _norm_M12_d,
							      _norm_M34_d,
							      _norm_CosTheta12_d,
							      _norm_CosTheta34_d,
							      _norm_phi_d,
							      0,
							      _norm_LS_d);
    } // end if computing cached values

    // do norm integral
    fptype normResult = NormEvents_4Body_Base::doNormIntegral_TD(
							  resolution,
							  tau,
							  xmixing,
							  ymixing,
							  dalitzId,
							  _totNumAccNormEvents,
							  _norm_SF_d,
							  _norm_LS_d);

    return normResult / _totNumAccNormEvents;
  }

} // end namespace GooFit
