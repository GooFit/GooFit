#include <numeric>
#include <vector>
#include <tuple>

#include <thrust/copy.h>

#include <mcbooster/GContainers.h>
#include <mcbooster/GTypes.h>

#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Weighted.h>
#include <goofit/MathUtils.h>

#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

namespace GooFit {

__host__ NormEvents_4Body_Weighted::NormEvents_4Body_Weighted(mcbooster::RealVector_h m12,
                                                              mcbooster::RealVector_h m34,
                                                              mcbooster::RealVector_h cos12,
                                                              mcbooster::RealVector_h cos34,
                                                              mcbooster::RealVector_h phi,
                                                              mcbooster::RealVector_h dtime,
                                                              mcbooster::RealVector_h sigma,
                                                              mcbooster::RealVector_h weights) {
    _norm_M12_h          = m12;
    _norm_M34_h          = m34;
    _norm_CosTheta12_h   = cos12;
    _norm_CosTheta34_h   = cos34;
    _norm_phi_h          = phi;
    _norm_dtime_h        = dtime;
    _norm_sigma_h        = sigma;
    _norm_weight_h       = weights;
    _totNumAccNormEvents = weights.size();
    _sumInitWeights      = thrust::reduce(thrust::host, weights.begin(), weights.end(), 0);

    _norm_SF_h = mcbooster::RealVector_h(0);
    _norm_LS_h = mcbooster::mc_host_vector<fpcomplex>(0);
}

__host__ fptype NormEvents_4Body_Weighted::computeNorm_TD(bool noCachedNormValuesToCompute,
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

    if(!noCachedNormValuesToCompute) {
        if(_norm_SF_h.size() != numSFCacheEntries) {
            _norm_SF_h.resize(numSFCacheEntries);
        }

        if(_norm_LS_h.size() != numLSCacheEntries) {
            _norm_LS_h.resize(numLSCacheEntries);
        }
    }

    // Copy cached SF values to the device.
    mcbooster::RealVector_d sf_batchResult_d = _norm_SF_h;
    // Copy cached lineshape values to the device.
    mcbooster::mc_device_vector<fpcomplex> ls_batchResult_d = _norm_LS_h;
    // Copy weights, dtime, and sigma to the device.
    mcbooster::RealVector_d batchWeights_d = _norm_weight_h;
    mcbooster::RealVector_d batchTime_d    = _norm_dtime_h;
    mcbooster::RealVector_d batchSigma_d   = _norm_sigma_h;

    if(!noCachedNormValuesToCompute) {
        // copy batch of norm events to device
        mcbooster::RealVector_d batchNormM12_d        = _norm_M12_h;
        mcbooster::RealVector_d batchNormM34_d        = _norm_M34_h;
        mcbooster::RealVector_d batchNormCosTheta12_d = _norm_CosTheta12_h;
        mcbooster::RealVector_d batchNormCosTheta34_d = _norm_CosTheta34_h;
        mcbooster::RealVector_d batchNormPhi_d        = _norm_phi_h;

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

        if(sfValuesUpdated) {
            _norm_SF_h = sf_batchResult_d;
        }

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

        if(lsValuesUpdated) {
            _norm_LS_h = ls_batchResult_d;
        }
    }

    return NormEvents_4Body_Base::doNormIntegral_TD(resolution,
                                                    tau,
                                                    xmixing,
                                                    ymixing,
                                                    dalitzId,
                                                    _totNumAccNormEvents,
                                                    sf_batchResult_d,
                                                    ls_batchResult_d,
                                                    batchTime_d,
                                                    batchSigma_d,
                                                    batchWeights_d);
}

} // namespace GooFit
