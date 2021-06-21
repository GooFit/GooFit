#pragma once

#include <goofit/PDFs/physics/Amp3BodyBase.h>
#include <goofit/PDFs/physics/Amp3Body_TD.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

class SpecialIncoherentIntegrator;
class SpecialIncoherentResonanceCalculator;

/**
Similar to `Amp3Body`, but the resonances
are added incoherently:

\f[
P(m^2_{12},m^2_{13};\vec\alpha) =
\sum\limits_i \left|\alpha_i B_i(m^2_{12},m^2_{13})\right|^2\epsilon(m^2_{12},m^2_{13})
\f]

The constructor is the same, but note that the `amp_imag` member of
GooFit::ResonancePdf is not used, so the \f$\alpha\f$ are in effect
interpreted as real numbers.

Very similar class to Amp3Body_TD, but without time dependence
(so no time resolution or mixing) and ignoring interference between
waves. This makes the code just different enough, the assumptions are
just enough changed, that it's not worth trying to modify or subclass
Amp3Body_TD to deal with both cases. So instead we have a separate
class with fairly similar structure.
**/

class Amp3Body_IS : public Amp3BodyBase {
  public:
    Amp3Body_IS(std::string n, Observable m12, Observable m13, EventNumber eventNumber, DecayInfo3 decay, GooPdf *eff);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // incoherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.
    __host__ auto normalize() -> fptype override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 3);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }

    __host__ void populateArrays() override;

  protected:
  private:
    DecayInfo3 decayInfo;
    Observable _m12;
    Observable _m13;
    fptype *dalitzNormRange;

    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<fpcomplex> *cachedResonances; // BW (and other) results for each event.
    double *integrals;                                  // Integrals of each BW resonance across the Daliz plot.

    bool *redoIntegral;
    mutable bool forceRedoIntegrals;
    fptype *cachedMasses;
    fptype *cachedWidths;
    int totalEventSize;
    int cacheToUse;
    PdfBase *efficiency;
    SpecialIncoherentIntegrator **integrators;
    SpecialIncoherentResonanceCalculator **calculators;
    int efficiencyFunction;
};

} // namespace GooFit
