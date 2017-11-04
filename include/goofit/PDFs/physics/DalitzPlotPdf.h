#pragma once

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

#include <thrust/complex.h>

namespace GooFit {

class SpecialResonanceIntegrator;
class SpecialResonanceCalculator;

class DalitzPlotPdf : public GooPdf {
  public:
    DalitzPlotPdf(
        std::string n, Variable *m12, Variable *m13, CountingVariable *eventNumber, DecayInfo *decay, GooPdf *eff);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalisation will get *really* confused and give wrong answers.

    __host__ fptype normalize() const override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 3);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }

    __host__ virtual void recursiveSetIndices();

  protected:
  private:
    DecayInfo *decayInfo;
    Variable *_m12;
    Variable *_m13;
    fptype *dalitzNormRange;

    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<thrust::complex<fptype>> *cachedWaves[16]; // Caches the BW values for each event.
    thrust::complex<fptype> ***integrals; // Caches the integrals of the BW waves for each combination of resonances.

    bool *redoIntegral;
    mutable bool forceRedoIntegrals;
    fptype *cachedMasses;
    fptype *cachedWidths;
    int totalEventSize;
    int cacheToUse;
    SpecialResonanceIntegrator ***integrators;
    SpecialResonanceCalculator **calculators;

    unsigned int efficiencyFunction;
};

class SpecialResonanceIntegrator
    : public thrust::unary_function<thrust::tuple<int, fptype *, int>, thrust::complex<fptype>> {
  public:
    // Class used to calculate integrals of terms BW_i * BW_j^*.
    SpecialResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj);
    void setDalitzIndex(unsigned int id) { dalitz_i = id; }
    void setResonanceIndex(unsigned int id) { resonance_i = id; }
    void setEfficiencyIndex(unsigned int id) { resonance_j = id; }
    __device__ thrust::complex<fptype> operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int dalitz_i;
    unsigned int resonance_i;
    unsigned int resonance_j;
    unsigned int parameters;
};

class SpecialResonanceCalculator
    : public thrust::unary_function<thrust::tuple<int, fptype *, int>, thrust::complex<fptype>> {
  public:
    // Used to create the cached BW values.
    SpecialResonanceCalculator(int pIdx, unsigned int res_idx);
    void setDalitzIndex(unsigned int id) { dalitz_i = id; }
    void setResonanceIndex(unsigned int id) { resonance_i = id; }
    __device__ thrust::complex<fptype> operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int dalitz_i;
    unsigned int resonance_i;
    unsigned int parameters;
};

} // namespace GooFit
