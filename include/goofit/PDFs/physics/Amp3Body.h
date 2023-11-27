#pragma once

#include <goofit/PDFs/physics/Amp3BodyBase.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

#include <mcbooster/GContainers.h>

#include <tuple>

#include <goofit/detail/Complex.h>

// Matrix
#include <Eigen/Core>
#include <Eigen/LU>

namespace GooFit {

class SpecialResonanceIntegrator;
class SpecialResonanceCalculator;
class DalitzPlotter;

/**
A time-independent description of the Dalitz plot
as a coherent sum of resonances:

\f[
P(m^2_{12},m^2_{13};\vec\alpha) = \left|\sum\limits_i
\alpha_i B_i(m^2_{12},m^2_{13})\right|^2
\epsilon(m^2_{12},m^2_{13})
\f]

where \f$\alpha_i\f$ is a complex coefficient, \f$B_i\f$ is a resonance
parametrization (see Goofit::ResonancePdf), and \f$\epsilon\f$ is a
real-valued efficiency function. The constructor takes the
squared-mass variables \f$m_{12}\f$ and \f$m_{13}\f$, an event index (this
is used in caching), a GooFit::DecayInfo3 object which contains a `vector`
of GooFit::ResonancePdf's as well as some global information like the
mother and daughter masses, and the efficiency function.
**/

class Amp3Body : public Amp3BodyBase {
  public:
    Amp3Body(std::string n,
             Observable m13,
             Observable m23,
             EventNumber eventNumber,
             DecayInfo3 decay,
             GooPdf *eff = nullptr);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    __host__ auto normalize() -> fptype override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 3, unsigned int offset = 0);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }

    __host__ void setGenerationOffset(int off) { generation_offset = off; }
    __host__ auto getGenerationOffset() -> int { return generation_offset; }
    __host__ auto GenerateSig(unsigned int numEvents, int seed = 0) -> std::
        tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h>;

    __host__ void populateArrays() override;

    /// Get the cached wave (device) vectors
    __host__ auto getCachedWaveNoCopy(size_t i) const -> const thrust::device_vector<fpcomplex> & {
        return *(cachedWaves[i]);
    }

    __host__ auto getCachedWave(size_t i) const -> const std::vector<std::complex<fptype>>;

    /// Sum up a cached wave
    __host__ auto sumCachedWave(size_t i) const -> fpcomplex;

    /// Get the decay info struct
    __host__ auto getDecayInfo() -> DecayInfo3 & { return decayInfo; }
    __host__ static void resetCacheCounter() { cacheCount = 0; }

    /// Calculate fit fractions (Cache should be pre-filled)
    __host__ auto fit_fractions(bool print = false) -> std::vector<std::vector<fptype>>;

    __host__ auto getFFTotalIntegral() const -> fptype { return totalFF_integral; }

    friend DalitzPlotter;

  protected:
    DecayInfo3 decayInfo;
    Observable _m13;
    Observable _m23;
    EventNumber _eventNumber;
    fptype *dalitzNormRange;

    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<fpcomplex> *cachedWaves[20]; // Caches the BW values for each event.
    fpcomplex ***integrals;    // Caches the integrals of the BW waves for each combination of resonances.
    fpcomplex ***integrals_ff; // Caches the integrals of the BW waves for each combination of resonances (for the fit
                               // fractions).

    mutable bool generation_no_norm{false};
    bool *redoIntegral;
    mutable bool forceRedoIntegrals;
    fptype *cachedMasses;
    fptype *cachedWidths;
    int totalEventSize;
    int eventOffset;
    int cacheToUse;
    static int cacheCount;
    int generation_offset{0};
    SpecialResonanceIntegrator ***integrators;
    SpecialResonanceIntegrator ***integrators_ff;
    SpecialResonanceCalculator **calculators;

    fptype totalFF_integral;

    unsigned int efficiencyFunction;
};

} // namespace GooFit
