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

std::vector<NormEvents_4Body_Base *>
NormEvents_4Body_HostCached::buildBatches(const std::vector<long> &normSeeds,
                                          unsigned int numNormEventsToGenPerBatch,
                                          const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses) {
    std::vector<NormEvents_4Body_Base *> ret(normSeeds.size());

    for(int n = 0; n < normSeeds.size(); n++) {
        ret[n] = new NormEvents_4Body_HostCached(normSeeds[n], numNormEventsToGenPerBatch, motherAndDaughterMasses);
    }

    return ret;
}

__host__ NormEvents_4Body_HostCached::NormEvents_4Body_HostCached(
    long normSeed,
    unsigned int numNormEventsToGenPerBatch,
    const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses) {
    // generate norm events on device then save to host
    mcbooster::RealVector_d temp_norm_M12_d(0);
    mcbooster::RealVector_d temp_norm_M34_d(0);
    mcbooster::RealVector_d temp_norm_CosTheta12_d(0);
    mcbooster::RealVector_d temp_norm_CosTheta34_d(0);
    mcbooster::RealVector_d temp_norm_phi_d(0);
    _totNumAccNormEvents = NormEvents_4Body_Base::generate4BodyNormEvents(normSeed,
                                                                          numNormEventsToGenPerBatch,
                                                                          motherAndDaughterMasses,
                                                                          temp_norm_M12_d,
                                                                          temp_norm_M34_d,
                                                                          temp_norm_CosTheta12_d,
                                                                          temp_norm_CosTheta34_d,
                                                                          temp_norm_phi_d);

    _norm_M12_h        = temp_norm_M12_d;
    _norm_M34_h        = temp_norm_M34_d;
    _norm_CosTheta12_h = temp_norm_CosTheta12_d;
    _norm_CosTheta34_h = temp_norm_CosTheta34_d;
    _norm_phi_h        = temp_norm_phi_d;

    _norm_SF_h = mcbooster::RealVector_h(0);
    _norm_LS_h = mcbooster::mc_host_vector<fpcomplex>(0);

    GOOFIT_INFO("Total # of accepted MC events used for normalization: {}", getNumAccNormEvents());
}

__host__ fptype NormEvents_4Body_HostCached::computeNorm_TD(bool noCachedNormValuesToCompute,
                                                            const MixingTimeResolution *const resolution,
                                                            fptype tau,
                                                            fptype xmixing,
                                                            fptype ymixing,
                                                            unsigned int dalitzId,
                                                            bool spinsCalculated,
                                                            const std::vector<bool> &lineshapeChanged,
                                                            const std::vector<unsigned int> &sfFunctionIndices,
                                                            const std::vector<unsigned int> &lsFunctionIndices) {
    unsigned int numSFCacheEntries = _totNumAccNormEvents * sfFunctionIndices.size();
    unsigned int numLSCacheEntries = _totNumAccNormEvents * lsFunctionIndices.size();

    // resize if needed
    if(!noCachedNormValuesToCompute) {
        if(_norm_SF_h.size() != numSFCacheEntries) {
            _norm_SF_h.resize(numSFCacheEntries);
        }

        if(_norm_LS_h.size() != numLSCacheEntries) {
            _norm_LS_h.resize(numLSCacheEntries);
        }
    }

    // do computations on device
    // copy previously cached sf values for batch to device
    mcbooster::RealVector_d sf_batchResult_d = _norm_SF_h;

    // copy previously cached ls values for batch to device
    mcbooster::mc_device_vector<fpcomplex> ls_batchResult_d = _norm_LS_h;

    if(!noCachedNormValuesToCompute) {
        // copy batch of norm events to device
        mcbooster::RealVector_d batchNormM12_d        = _norm_M12_h;
        mcbooster::RealVector_d batchNormM34_d        = _norm_M34_h;
        mcbooster::RealVector_d batchNormCosTheta12_d = _norm_CosTheta12_h;
        mcbooster::RealVector_d batchNormCosTheta34_d = _norm_CosTheta34_h;
        mcbooster::RealVector_d batchNormPhi_d        = _norm_phi_h;

        // recompute cached sf values for batch (if needed)
        bool sfValuesUpdated = NormEvents_4Body_Base::computeCachedSFValuesForBatch_TD(spinsCalculated,
                                                                                       dalitzId,
                                                                                       sfFunctionIndices,
                                                                                       _totNumAccNormEvents,
                                                                                       batchNormM12_d,
                                                                                       batchNormM34_d,
                                                                                       batchNormCosTheta12_d,
                                                                                       batchNormCosTheta34_d,
                                                                                       batchNormPhi_d,
                                                                                       0,
                                                                                       sf_batchResult_d);

        // copy sf values for batch back to host if any of the values changed
        if(sfValuesUpdated) // TODO only copy changed values
        {
            _norm_SF_h = sf_batchResult_d;
        }

        // recompute cached ls values for batch (if needed)
        bool lsValuesUpdated = NormEvents_4Body_Base::computeCachedLSValuesForBatch_TD(lineshapeChanged,
                                                                                       dalitzId,
                                                                                       lsFunctionIndices,
                                                                                       _totNumAccNormEvents,
                                                                                       batchNormM12_d,
                                                                                       batchNormM34_d,
                                                                                       batchNormCosTheta12_d,
                                                                                       batchNormCosTheta34_d,
                                                                                       batchNormPhi_d,
                                                                                       0,
                                                                                       ls_batchResult_d);

        // copy ls values back to host if any of the values changed
        if(lsValuesUpdated) // TODO only copy changed values
        {
            _norm_LS_h = ls_batchResult_d;
        }
    } // end if computing cached values

    // do norm integral
    return NormEvents_4Body_Base::doNormIntegral_TD(
        resolution, tau, xmixing, ymixing, dalitzId, _totNumAccNormEvents, sf_batchResult_d, ls_batchResult_d);
}

} // end namespace GooFit
