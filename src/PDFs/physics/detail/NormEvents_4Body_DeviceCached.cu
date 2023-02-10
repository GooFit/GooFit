#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_DeviceCached.h>

namespace GooFit {

std::vector<NormEvents_4Body_Base *>
NormEvents_4Body_DeviceCached::buildBatches(const std::vector<long> &normSeeds,
                                            unsigned int numNormEventsToGenPerBatch,
                                            const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses) {
    std::vector<NormEvents_4Body_Base *> ret(normSeeds.size());

    for(int n = 0; n < normSeeds.size(); n++) {
        ret[n] = new NormEvents_4Body_DeviceCached(motherAndDaughterMasses, normSeeds[n], numNormEventsToGenPerBatch);
    }

    return ret;
}

__host__ NormEvents_4Body_DeviceCached::NormEvents_4Body_DeviceCached(
    const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses, long normSeed, unsigned int numNormEventsToGen) {
    _norm_M12_d        = mcbooster::RealVector_d(0);
    _norm_M34_d        = mcbooster::RealVector_d(0);
    _norm_CosTheta12_d = mcbooster::RealVector_d(0);
    _norm_CosTheta34_d = mcbooster::RealVector_d(0);
    _norm_phi_d        = mcbooster::RealVector_d(0);

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
    // Events are unweighted.
    _sumInitWeights = _totNumAccNormEvents;

    GOOFIT_INFO("# of accepted MC events used for normalization: {}", getNumAccNormEvents());
}

__host__ fptype NormEvents_4Body_DeviceCached::computeNorm_TD(bool noCachedNormValuesToCompute,
                                                              const MixingTimeResolution *const resolution,
                                                              fptype tau,
                                                              fptype xmixing,
                                                              fptype ymixing,
                                                              unsigned int dalitzId,
                                                              bool spinsCalculated,
                                                              const std::vector<bool> &lineshapeChanged,
                                                              const std::vector<unsigned int> &sfFunctionIndices,
                                                              const std::vector<unsigned int> &lsFunctionIndices) {
    if(!noCachedNormValuesToCompute) {
        unsigned int numSFCacheEntries = _totNumAccNormEvents * sfFunctionIndices.size();
        if(_norm_SF_d.size() != numSFCacheEntries) {
            _norm_SF_d.resize(numSFCacheEntries);
        }

        // compute cached sf values (if needed)
        NormEvents_4Body_Base::computeCachedSFValuesForBatch_TD(spinsCalculated,
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
        if(_norm_LS_d.size() != numLSCacheEntries) {
            _norm_LS_d.resize(numLSCacheEntries);
        }

        // recompute cached ls values for batch (if needed)
        NormEvents_4Body_Base::computeCachedLSValuesForBatch_TD(lineshapeChanged,
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
        resolution, tau, xmixing, ymixing, dalitzId, _totNumAccNormEvents, _norm_SF_d, _norm_LS_d);

    return normResult;
}

} // end namespace GooFit
