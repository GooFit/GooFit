/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include <goofit/PDFs/physics/Amp4BodyBase.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/NormSpinCalculator_TD.h>

#include <mcbooster/GContainers.h>

#include <thrust/remove.h>

#include <tuple>

namespace GooFit {

class LSCalculator_TD;
class AmpCalc_TD;
class SFCalculator_TD;
class NormIntegrator_TD;
class Lineshape;

class Amp4Body_TD final : public Amp4BodyBase {
  public:
    Amp4Body_TD(std::string n,
                std::vector<Observable> observables,
                DecayInfo4t decay,
                MixingTimeResolution *r,
                GooPdf *eff,
                Observable *mistag,
		long normSeed,
                unsigned int MCeventsNorm);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    __host__ fptype normalize() override;

    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 8);
    __host__ void setForceIntegrals(bool f = true) { _forceRedoIntegrals = f; }
    __host__ int getNumAccNormEvents() const { return _nAcc_Norm_Events; } 
    __host__ void setGenerationOffset(int off) { generation_offset = off; }
    __host__ void setMaxWeight(fptype wmax) { maxWeight = wmax; }

    __host__ std::
        tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
        GenerateSig(unsigned int numEvents, int seed = 0);

    __host__ void populateArrays() override;

    void printAmpMappings() const;

    void printSelectedLineshapes(const std::vector<unsigned int>& lsIndices) const;

    void printSelectedSFs(const std::vector<unsigned int>& sfIndices) const;

  protected:
  private:
    __host__ void computeCachedNormValues(const std::vector<bool>& lineshapeChanged);

    __host__ void computeCachedValues(const std::vector<bool>& lineshapeChanged, const std::vector<bool>& amplitudeComponentChanged);

    __host__ std::vector<bool> areLineshapesChanged() const;

    __host__ std::vector<bool> areAmplitudeComponentsChanged() const;

    __host__ void generateAndSetNormEvents();

    std::vector<SpinFactor *> _SpinFactors;
    std::vector<Lineshape *> _LineShapes;
    std::vector<AmpCalc_TD *> _AmpCalcs;
    NormIntegrator_TD *_Integrator;
    std::vector<SFCalculator_TD *> _sfcalculators;
    std::vector<LSCalculator_TD *> _lscalculators;

    unsigned int efficiencyFunction;

    // store normalization events
    mcbooster::RealVector_d _norm_M12;
    mcbooster::RealVector_d _norm_M34;
    mcbooster::RealVector_d _norm_CosTheta12;
    mcbooster::RealVector_d _norm_CosTheta34;
    mcbooster::RealVector_d _norm_phi;
    // store spin and lineshape values for normalization
    mutable mcbooster::RealVector_d _norm_SF;
    mutable mcbooster::mc_device_vector<fpcomplex> _norm_LS;

    DecayInfo4t _decayInfo;
    MixingTimeResolution *resolution;
    int _nAcc_Norm_Events;
    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<fpcomplex> *_cachedResSF{nullptr}; // Caches the BW values and Spins for each event.
    thrust::device_vector<fpcomplex> *_cachedAMPs{nullptr};  // cache Amplitude values for each event.
    mutable bool _generation_no_norm{false};
    mutable bool _SpinsCalculated{false};
    mutable bool _forceRedoIntegrals{true};
    int _totalEventSize;
    int cacheToUse{0};
    unsigned int generation_offset{0};
    double maxWeight{0};
    unsigned int _NUM_AMPLITUDES;
    const long _NORM_SEED;
    const int _NUM_NORM_EVENTS_TO_GEN;
};

} // namespace GooFit
