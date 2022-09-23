#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_DeviceCached.h>
#include<iostream>
#include <cstdlib>


namespace GooFit {

std::vector<NormEvents_4Body_Base *>
NormEvents_4Body_DeviceCached::buildBatches(const std::vector<long> &normSeeds,
                                            unsigned int numNormEventsToGenPerBatch,
                                            const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses,
                                            bool with_acceptance) {
    std::vector<NormEvents_4Body_Base *> ret(normSeeds.size());

    for(int n = 0; n < normSeeds.size(); n++) {
        NormEvents_4Body_DeviceCached* device_batch = new NormEvents_4Body_DeviceCached(motherAndDaughterMasses, normSeeds[n], numNormEventsToGenPerBatch);
        if(true){
            //std::vector<mcbooster::RealVector_h> norm_phsp = device_batch->get_norm_phsp();
            //mcbooster::RealVector_h test_norm, test_dtime, test_dtime_weights;
            //device_batch->set_norm_info(test_dtime, test_dtime_weights, test_norm );

            device_batch->write_norm_phsp("test_norm_events.txt");
            //execute python script here
            std::system("python assign_acceptance_weights.py");
            //check that python script created weights file and it is not empty!
            device_batch->read_norm_info("norm_weights.txt");
            std::vector<mcbooster::RealVector_h> norm_info = device_batch->get_norm_info();
            for(int i = 0; i < 10;i++){
                printf("dtime: %.7g, dtime_weight: %.7g, eff_weight: %.7g\n",norm_info[0][i], norm_info[1][i], norm_info[2][i]);
            }
        }
        
        //for(int i = 0; i < 10;i++){
        //    printf("m12: %.7g ",norm_phsp[0][i]);
        //}
        //printf("\n");

        ret[n] = device_batch;
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
                                                              const std::vector<unsigned int> &lsFunctionIndices,
                                                              unsigned int CacheIdx) {
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
        resolution, tau, xmixing, ymixing, dalitzId, _totNumAccNormEvents, _norm_SF_d, _norm_LS_d, CacheIdx);

    return normResult;
}

} // end namespace GooFit
