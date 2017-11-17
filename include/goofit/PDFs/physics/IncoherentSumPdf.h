#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/TddpPdf.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

// Very similar class to TddpPdf, but without time dependence
// (so no time resolution or mixing) and ignoring interference between
// waves. This makes the code just different enough, the assumptions are
// just enough changed, that it's not worth trying to modify or subclass
// TddpPdf to deal with both cases. So instead we have a separate
// class with fairly similar structure.

class SpecialIncoherentIntegrator;
class SpecialIncoherentResonanceCalculator;

class IncoherentSumPdf : public GooPdf {
  public:
    IncoherentSumPdf(
        std::string n, Observable m12, Observable m13, EventNumber eventNumber, DecayInfo3 decay, GooPdf *eff);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // incoherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalisation will get *really* confused and give wrong answers.
    __host__ fptype normalize() const override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 3);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }

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
};

class SpecialIncoherentIntegrator : public thrust::unary_function<thrust::tuple<int, fptype *>, fptype> {
  public:
    SpecialIncoherentIntegrator(int pIdx, unsigned int ri);
    __device__ fptype operator()(thrust::tuple<int, fptype *> t) const;

  private:
    unsigned int resonance_i;
    unsigned int parameters;
};

class SpecialIncoherentResonanceCalculator
    : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fpcomplex> {
  public:
    SpecialIncoherentResonanceCalculator(int pIdx, unsigned int res_idx);
    __device__ fpcomplex operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int resonance_i;
    unsigned int parameters;
};

} // namespace GooFit
