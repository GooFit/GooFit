#include <numeric>
#include <vector>
#include <tuple>

#include <thrust/copy.h>

#include <mcbooster/GContainers.h>
#include <mcbooster/GTypes.h>

#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_WeightedDevice.h>
#include <goofit/MathUtils.h>

#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

namespace GooFit {

__host__ NormEvents_4Body_WeightedDevice::NormEvents_4Body_WeightedDevice(mcbooster::RealVector_d m12,
                                                                          mcbooster::RealVector_d m34,
                                                                          mcbooster::RealVector_d cos12,
                                                                          mcbooster::RealVector_d cos34,
                                                                          mcbooster::RealVector_d phi,
                                                                          mcbooster::RealVector_d dtime,
                                                                          mcbooster::RealVector_d sigma,
                                                                          mcbooster::RealVector_d weights) {
    _norm_M12_d          = m12;
    _norm_M34_d          = m34;
    _norm_CosTheta12_d   = cos12;
    _norm_CosTheta34_d   = cos34;
    _norm_phi_d          = phi;
    _norm_dtime_d        = dtime;
    _norm_sigma_d        = sigma;
    _norm_weight_d       = weights;
    _totNumAccNormEvents = weights.size();
    _sumInitWeights      = thrust::reduce(thrust::device, weights.begin(), weights.end(), 0);

    _norm_SF_d = mcbooster::RealVector_d(0);
    _norm_LS_d = mcbooster::mc_device_vector<fpcomplex>(0);
}

__host__ fptype NormEvents_4Body_WeightedDevice::computeNorm_TD(bool noCachedNormValuesToCompute,
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
        if(_norm_SF_d.size() != numSFCacheEntries) {
            _norm_SF_d.resize(numSFCacheEntries);
        }

        if(_norm_LS_d.size() != numLSCacheEntries) {
            _norm_LS_d.resize(numLSCacheEntries);
        }
    }

    if(!noCachedNormValuesToCompute) {
        bool sfValuesUpdated = NormEvents_4Body_Base::computeCachedSFValuesForBatch_TD(spinsCalculated,
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

        bool lsValuesUpdated = NormEvents_4Body_Base::computeCachedLSValuesForBatch_TD(lineshapeChanged,
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
    }

    return NormEvents_4Body_Base::doNormIntegral_TD(resolution,
                                                    tau,
                                                    xmixing,
                                                    ymixing,
                                                    dalitzId,
                                                    _totNumAccNormEvents,
                                                    _norm_SF_d,
                                                    _norm_LS_d,
                                                    _norm_dtime_d,
                                                    _norm_sigma_d,
                                                    _norm_weight_d);
}

} // namespace GooFit
