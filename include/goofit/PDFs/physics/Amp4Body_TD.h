/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficiently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include <memory>

#include <goofit/PDFs/physics/Amp4BodyBase.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/NormSpinCalculator_TD.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>

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
    Amp4Body_TD()                                     = delete;
    Amp4Body_TD(const Amp4Body_TD &copyMe)            = delete;
    Amp4Body_TD(Amp4Body_TD &&moveMe)                 = delete;
    Amp4Body_TD &operator=(const Amp4Body_TD &copyMe) = delete;
    Amp4Body_TD &operator=(Amp4Body_TD &&moveMe)      = delete;
    virtual ~Amp4Body_TD() override                   = default;

    // Build Amp4Body_TD where the MC events used for normalization are stored on the host side
    // and where normalization computations are done in batches.
    // This is much slower than the other case (where the events used for normalization are stored on the device side)
    // but allows you to use a larger # of events for normalization than would otherwise would be possible.
    __host__ Amp4Body_TD(std::string n,
                         std::vector<Observable> observables,
                         DecayInfo4t decay,
                         MixingTimeResolution *Tres,
                         GooPdf *efficiency,
                         Observable *mistag,
                         const std::vector<long> &normSeeds,
                         unsigned int numNormEventsToGenPerBatch);

    // Build Amp4Body_TD where the MC events used for normalization are stored on the device side.
    // This is much faster than the other case (where the events used for normalization are stored on the host side)
    // but places a lower limit on the maximum # of events that can be used for normalization.
    __host__ Amp4Body_TD(std::string n,
                         std::vector<Observable> observables,
                         DecayInfo4t decay,
                         MixingTimeResolution *Tres,
                         GooPdf *efficiency,
                         Observable *mistag,
                         long normSeed,
                         unsigned int numNormEventsToGen);

    // Does common initialization
    __host__ Amp4Body_TD(std::string n,
                         std::vector<Observable> observables,
                         DecayInfo4t decay,
                         MixingTimeResolution *Tres,
                         GooPdf *efficiency,
                         Observable *mistag,
                         const std::vector<NormEvents_4Body_Base *> &normEvents);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    __host__ auto normalize() -> fptype override;

    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 8);

    __host__ void setForceIntegrals(bool f = true) { _forceRedoIntegrals = f; }

    __host__ int getNumAccNormEvents() const;

    __host__ fptype getSumInitNormEventWeights() const;

    __host__ void setGenerationOffset(int off) { generation_offset = off; }

    __host__ void setMaxWeight(fptype wmax) { maxWeight = wmax; }

    __host__ auto GenerateSig(unsigned int numEvents, int seed = 0) -> std::
        tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>;

    __host__ void populateArrays() override;

    void printAmpMappings() const;

    void printSelectedLineshapes(const std::vector<unsigned int> &lsIndices) const;

    void printSelectedSFs(const std::vector<unsigned int> &sfIndices) const;

  protected:
  private:
    __host__ void computeCachedValues(const std::vector<bool> &lineshapeChanged,
                                      const std::vector<bool> &amplitudeComponentChanged);

    __host__ std::vector<bool> areLineshapesChanged() const;

    __host__ std::vector<bool> areAmplitudeComponentsChanged() const;

    __host__ std::vector<unsigned int> getSFFunctionIndices() const;

    __host__ std::vector<unsigned int> getLSFunctionIndices() const;

    std::vector<SpinFactor *> _SpinFactors;
    std::vector<Lineshape *> _LineShapes;
    std::vector<AmpCalc_TD *> _AmpCalcs;
    std::vector<SFCalculator_TD *> _sfcalculators;
    std::vector<LSCalculator_TD *> _lscalculators;
    unsigned int efficiencyFunction;
    std::vector<std::unique_ptr<NormEvents_4Body_Base>> _normEvents;
    const DecayInfo4t _DECAY_INFO;
    MixingTimeResolution *_resolution;
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
};

} // namespace GooFit
