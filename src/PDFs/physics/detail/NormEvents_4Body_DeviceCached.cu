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


  __host__ NormEvents_4Body_Batch NormEvents_4Body_DeviceCached::getBatch(unsigned int batchNum) const
  {
    if (batchNum >= _numBatches)
    {
      throw GooFit::GeneralError("Invalid batch number for NormEvents_4Body_DeviceCached object.");
    }

    return NormEvents_4Body_Batch(_norm_M12_d, 
				  _norm_M34_d,
				  _norm_CosTheta12_d, 
				  _norm_CosTheta34_d,
				  _norm_phi_d,
				  getNumAccNormEvents(),
				  0);
  }

} // end namespace GooFit
