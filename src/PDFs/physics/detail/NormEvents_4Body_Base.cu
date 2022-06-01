#include <vector>

#include <thrust/transform_reduce.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/remove.h>

#include <mcbooster/Generate.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/EvaluateArray.h>

#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/PDFs/physics/detail/NormLSCalculator_TD.h>
#include <goofit/PDFs/physics/detail/NormSpinCalculator_TD.h>
#include <goofit/PDFs/physics/detail/FourDblTupleAdd.h>
#include <goofit/PDFs/physics/detail/NormIntegrator_TD.h>

namespace GooFit {

__host__ unsigned int
NormEvents_4Body_Base::generate4BodyNormEvents(long normSeed,
                                               unsigned int numNormEventsToGen,
                                               const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses,
                                               mcbooster::RealVector_d &norm_M12,
                                               mcbooster::RealVector_d &norm_M34,
                                               mcbooster::RealVector_d &norm_CosTheta12,
                                               mcbooster::RealVector_d &norm_CosTheta34,
                                               mcbooster::RealVector_d &norm_phi) {
    std::vector<mcbooster::GReal_t> masses(motherAndDaughterMasses.cbegin() + 1, motherAndDaughterMasses.cend());
    mcbooster::PhaseSpace phsp(motherAndDaughterMasses[0], masses, numNormEventsToGen, normSeed);
    phsp.Generate(mcbooster::Vector4R(motherAndDaughterMasses[0], 0.0, 0.0, 0.0));
    phsp.Unweight();

    auto nAcc                     = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1                       = phsp.GetDaughters(0);
    auto d2                       = phsp.GetDaughters(1);
    auto d3                       = phsp.GetDaughters(2);
    auto d4                       = phsp.GetDaughters(3);

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    auto zip_end   = zip_begin + d1.size();
    auto new_end   = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    d1.erase(thrust::get<0>(new_end.get_iterator_tuple()), d1.end());
    d2.erase(thrust::get<1>(new_end.get_iterator_tuple()), d2.end());
    d3.erase(thrust::get<2>(new_end.get_iterator_tuple()), d3.end());
    d4.erase(thrust::get<3>(new_end.get_iterator_tuple()), d4.end());

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    norm_M12.resize(nAcc);
    norm_M34.resize(nAcc);
    norm_CosTheta12.resize(nAcc);
    norm_CosTheta34.resize(nAcc);
    norm_phi.resize(nAcc);

    mcbooster::VariableSet_d VarSet(5);
    VarSet[0] = &norm_M12;
    VarSet[1] = &norm_M34;
    VarSet[2] = &norm_CosTheta12;
    VarSet[3] = &norm_CosTheta34;
    VarSet[4] = &norm_phi;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet);

    phsp.FreeResources();

    return nAcc;
}

__host__ bool
NormEvents_4Body_Base::computeCachedLSValuesForBatch_TD(const std::vector<bool> &lineshapeChanged,
                                                        unsigned int dalitzId,
                                                        const std::vector<unsigned int> &lsFunctionIndices,
                                                        unsigned int numAccThisBatch,
                                                        const mcbooster::RealVector_d &batchNormM12_d,
                                                        const mcbooster::RealVector_d &batchNormM34_d,
                                                        const mcbooster::RealVector_d &batchNormCosTheta12_d,
                                                        const mcbooster::RealVector_d &batchNormCosTheta34_d,
                                                        const mcbooster::RealVector_d &batchNormPhi_d,
                                                        unsigned int resultOffset,
                                                        mcbooster::mc_device_vector<fpcomplex> &ls_batchResult_d) {
    bool lsValuesUpdated = false;

    // lineshape value calculation for the normalization, also recalculated every time parameter change
    for(int i = 0; i < lsFunctionIndices.size(); ++i) {
        if(!lineshapeChanged[i]) {
            continue;
        }

        lsValuesUpdated = true;

        NormLSCalculator_TD ns;
        ns.setDalitzId(dalitzId);
        ns.setResonanceId(lsFunctionIndices[i]);

        thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(batchNormM12_d.cbegin(),
                                                                       batchNormM34_d.cbegin(),
                                                                       batchNormCosTheta12_d.cbegin(),
                                                                       batchNormCosTheta34_d.cbegin(),
                                                                       batchNormPhi_d.cbegin())),
                          thrust::make_zip_iterator(thrust::make_tuple(batchNormM12_d.cend(),
                                                                       batchNormM34_d.cend(),
                                                                       batchNormCosTheta12_d.cend(),
                                                                       batchNormCosTheta34_d.cend(),
                                                                       batchNormPhi_d.cend())),
                          (ls_batchResult_d.begin() + resultOffset + (i * numAccThisBatch)),
                          ns);
    }

    return lsValuesUpdated;
}

__host__ bool
NormEvents_4Body_Base::computeCachedSFValuesForBatch_TD(bool spinsCalculated,
                                                        unsigned int dalitzId,
                                                        const std::vector<unsigned int> &sfFunctionIndices,
                                                        unsigned int numAccThisBatch,
                                                        const mcbooster::RealVector_d &batchNormM12_d,
                                                        const mcbooster::RealVector_d &batchNormM34_d,
                                                        const mcbooster::RealVector_d &batchNormCosTheta12_d,
                                                        const mcbooster::RealVector_d &batchNormCosTheta34_d,
                                                        const mcbooster::RealVector_d &batchNormPhi_d,
                                                        unsigned int resultOffset,
                                                        mcbooster::RealVector_d &sf_batchResult_d) {
    bool sfValuesUpdated = false;

    // Calculate spinfactors only once for normalization events and real events
    // strided_range is a template implemented in DalitsPlotHelpers.hh
    // it basically goes through the array by increasing the pointer by a certain amount instead of just one step.
    if(!spinsCalculated) {
        for(int i = 0; i < sfFunctionIndices.size(); ++i) {
            sfValuesUpdated = true;

            NormSpinCalculator_TD nsc;
            nsc.setDalitzId(dalitzId);
            nsc.setSpinFactorId(sfFunctionIndices[i]);

            thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(batchNormM12_d.cbegin(),
                                                                           batchNormM34_d.cbegin(),
                                                                           batchNormCosTheta12_d.cbegin(),
                                                                           batchNormCosTheta34_d.cbegin(),
                                                                           batchNormPhi_d.cbegin())),
                              thrust::make_zip_iterator(thrust::make_tuple(batchNormM12_d.cend(),
                                                                           batchNormM34_d.cend(),
                                                                           batchNormCosTheta12_d.cend(),
                                                                           batchNormCosTheta34_d.cend(),
                                                                           batchNormPhi_d.cend())),
                              (sf_batchResult_d.begin() + resultOffset + (i * numAccThisBatch)),
                              nsc);
        }
    }

    return sfValuesUpdated;
}

__host__ fptype NormEvents_4Body_Base::doNormIntegral_TD(
    const MixingTimeResolution *const resolution,
    fptype tau,
    fptype xmixing,
    fptype ymixing,
    unsigned int dalitzId,
    unsigned int numAccThisBatch,
    mcbooster::RealVector_d
        &batchSF_d, // not modified, can't seem to make const type work with thrust::transform_reduce
    mcbooster::mc_device_vector<fpcomplex>
        &batchLS_d) // not modified, can't seem to make const type work with thrust::transform_reduce
{
    thrust::constant_iterator<fptype *> normSFaddress(thrust::raw_pointer_cast(batchSF_d.data()));
    thrust::constant_iterator<fpcomplex *> normLSaddress(thrust::raw_pointer_cast(batchLS_d.data()));
    thrust::constant_iterator<int> NumNormEvents(numAccThisBatch);
    thrust::counting_iterator<int> eventIndex(0);

    thrust::tuple<fptype, fptype, fptype, fptype> dummy(0, 0, 0, 0);
    FourDblTupleAdd MyFourDoubleTupleAdditionFunctor;
    thrust::tuple<fptype, fptype, fptype, fptype> sumIntegral;

    NormIntegrator_TD integrator;
    integrator.setDalitzId(dalitzId);

    sumIntegral = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, NumNormEvents, normSFaddress, normLSaddress)),
        thrust::make_zip_iterator(
            thrust::make_tuple(eventIndex + numAccThisBatch, NumNormEvents, normSFaddress, normLSaddress)),
        integrator,
        dummy,
        MyFourDoubleTupleAdditionFunctor);

    // GOOFIT_TRACE("sumIntegral={}", sumIntegral);

    // printf("normalize A2/#evts , B2/#evts: %.5g, %.5g\n",thrust::get<0>(sumIntegral)/_nAcc_Norm_Events,

    return resolution->normalization(thrust::get<0>(sumIntegral),
                                     thrust::get<1>(sumIntegral),
                                     thrust::get<2>(sumIntegral),
                                     thrust::get<3>(sumIntegral),
                                     tau,
                                     xmixing,
                                     ymixing);
}

} // end namespace GooFit
