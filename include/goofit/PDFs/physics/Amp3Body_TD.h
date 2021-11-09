#pragma once

#include <goofit/PDFs/physics/Amp3BodyBase.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

namespace GooFit {

class SpecialDalitzIntegrator;
class SpecialWaveCalculator;

/**
If the Gaussian is a potato, this is a five-course
banquet dinner involving entire roasted animals stuffed with other
animals, large dance troupes performing between the courses, an
orchestra playing in the background, and lengthy speeches. There
will not be a vegetarian option. Without going too deeply into the
physics, the function models a decay, eg \f$D^0\to\pi\pi\pi^0\f$, that
can happen either directly or through a mixed path
\f$D^0\to \overline{D^0}\to\pi\pi\pi^0\f$. (Although developed for the
\f$\pi\pi\pi^0\f$ case, it should be useful for any decay where the
final state is its own anti-state.) The probability of the mixing
path depends on the decay time, and quantum-mechanically interferes
with the direct path. Consequently the full Time-Dependent
Dalitz-Plot (Tddp) amplitude is (suppressing the dependence on
squared masses, for clarity):

\f[
\label{eq:fullmix}
\begin{align}
P(m^2_{12}, m^2_{13}, t, \sigma_t;x,y,\tau,\vec\alpha) &=&
e^{-t/\tau}\Big(|A+B|^2\cosh(yt/\tau)\\
&& + |A-B|^2\cos(xt/\tau)\\
&& - 2\Re(AB^*)\sinh(yt/\tau)\\
&& - 2\Im(AB^*)\sin(xt/\tau)\Big)
\end{align}
\f]

where (notice the
reversed masses in the \f$B\f$ calculation)

\f[
\begin{align}
A &=& \sum\limits_i \alpha_iB_i(m^2_{12}, m^2_{13}) \\
B &=& \sum\limits_i \alpha_iB_i(m^2_{13}, m^2_{12}),
\end{align}
\f]

*convolved with* a time-resolution function and *multiplied by* an
efficiency. The implementation involves a large amount of caching of
the intermediate \f$B_i\f$ values, because these are expected to change
slowly relative to the coefficients \f$\alpha\f$ (in many cases, not at
all, since masses and widths are often held constant) and are
relatively expensive to calculate.

The constructor takes the measured decay time \f$t\f$, error on decay
time \f$\sigma_t\f$, squared masses \f$m^2_{12}\f$ and \f$m^2_{13}\f$, event
number, decay information (the same class as in `Amp3Body`; it
also holds the mixing parameters \f$x\f$ and \f$y\f$ and lifetime \f$\tau\f$),
time-resolution function, efficiency, and optionally a mistag
fraction. A variant constructor takes, instead of a single
time-resolution function, a `vector` of functions and an additional
observable \f$m_{D^0}\f$; in this case the resolution function used
depends on which bin of \f$m_{D^0}\f$ the event is in, and the number of
bins is taken as equal to the number of resolution functions
supplied.

It is not suggested to try to use this thing from scratch. Start
with a working example and modify it gradually.
**/

class Amp3Body_TD : public Amp3BodyBase {
  public:
    Amp3Body_TD(std::string n,
                Observable _dtime,
                Observable _sigmat,
                Observable m12,
                Observable m13,
                EventNumber eventNumber,
                DecayInfo3t decay,
                MixingTimeResolution *r,
                GooPdf *eff,
                Observable *mistag   = nullptr,
                Observable *charmtag = nullptr);
    Amp3Body_TD(std::string n,
                Observable _dtime,
                Observable _sigmat,
                Observable m12,
                Observable m13,
                EventNumber eventNumber,
                DecayInfo3t decay,
                std::vector<MixingTimeResolution *> &r,
                GooPdf *eff,
                Observable md0,
                Observable *mistag   = nullptr,
                Observable *charmtag = nullptr);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    // The mistag variable is the probability that an event has a mother particle
    // that was correctly reconstructed but wrongly tagged. Consider an analysis
    // with three components: Signal, mistagged signal, and background. We want to
    // have the PDF be a sum, thus:
    // P = p_s S(m+, m-) + p_m(l_f S(m+, m-) + (1 - l_f)S(m-, m+)) + p_B B(m+, m-)
    // where p_s, p_m, p_B are the respective probabilities that this event
    // are signal, mistagged, or background, and l_f is the "lucky fraction",
    // that fraction of the mistagged signal which, by chance, got assigned
    // the correct charge. ('Mistagged' means that the wrong track was used
    // to determine charge, but about 50% of random tracks will have the same
    // charge as the right track did.) Clearly the above can be simplified (using
    // S and S' to indicate non-flipped and flipped versions of the signal):
    // P = (p_s + p_m*l_f) S + p_m*(1-l_f) S' + p_B B
    //   = a(bS + (1-b)S') + p_B B.
    // where
    // a = p_s + p_m
    // b = (p_s + p_m*l_f) / (p_s + p_m)
    // or in other words, b is the fraction of signal + mistag that has the right
    // charge. It's up to the user to create this variable. The default is to take
    // b as 1, if 'mistag' is not supplied.
    // Note that normalization is not affected because the integrals of S and S'
    // are identical and the weights sum to one, and efficiency is not affected
    // because it depends on the momenta of the daughter tracks, which are not
    // affected by making the wrong charge assignment to the mother.

    __host__ auto normalize() -> fptype override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 5, unsigned int offset = 0);
    __host__ void setD0Fraction(fptype d0fraction);
    __host__ fptype getD0Fraction();
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }
    __host__ static void resetCacheCounter() { cacheCount = 0; }

    /// Get the decay info struct
    __host__ DecayInfo3t &getDecayInfo() { return decayInfo; }
    /// Get the cached wave (device) vectors
    __host__ const thrust::device_vector<WaveHolder_s> &getCachedWaveNoCopy(size_t i) const {
        return *(cachedWaves[i]);
    }
    /// Calculate fit fractions (Cache should be pre-filled)
    __host__ std::vector<std::vector<fptype>> getFractions();

    __host__ void populateArrays() override;

  protected:
  private:
    DecayInfo3t decayInfo;
    Observable _m12;
    Observable _m13;
    MixingTimeResolution *resolution;
    GooPdf *_efficiency;
    Observable _mistag;
    Observable _charmtag;
    fptype *dalitzNormRange{nullptr};

    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<WaveHolder_s> *cachedWaves[16]; // Caches the BW values for each event.
    ThreeComplex ***integrals{nullptr}; // Caches the integrals of the BW waves for each combination of resonances.

    bool *redoIntegral;
    mutable bool forceRedoIntegrals{true};
    fptype *cachedMasses;
    fptype *cachedWidths;

    unsigned int resolutionFunction;
    unsigned int efficiencyFunction;

    int totalEventSize;
    int eventOffset;
    int cacheToUse{0};
    static int cacheCount;
    SpecialDalitzIntegrator ***integrators{nullptr};
    SpecialWaveCalculator **calculators{nullptr};

    fptype _D0Fraction;
};

} // namespace GooFit
