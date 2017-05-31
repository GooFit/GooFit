#pragma once

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"
#include "goofit/PDFs/physics/MixingTimeResolution_Aux.h"

namespace GooFit {

// thrust::tuple can't go down the read-only cache pipeline, so we are creating a structure for this.
typedef struct {
    // note: combining these into a single transaction (double2) should improve memory performance
    fptype ai_real;
    fptype ai_imag;
    fptype bi_real;
    fptype bi_imag;
} WaveHolder_s;

typedef thrust::tuple<fptype, fptype, fptype, fptype> WaveHolder;
typedef thrust::tuple<fptype, fptype, fptype, fptype, fptype, fptype> ThreeComplex;

class SpecialDalitzIntegrator;
class SpecialWaveCalculator;

class TddpPdf : public GooPdf {
  public:
    TddpPdf(std::string n,
            Variable *_dtime,
            Variable *_sigmat,
            Variable *m12,
            Variable *m13,
            CountingVariable *eventNumber,
            DecayInfo *decay,
            MixingTimeResolution *r,
            GooPdf *eff,
            Variable *mistag = nullptr);
    TddpPdf(std::string n,
            Variable *_dtime,
            Variable *_sigmat,
            Variable *m12,
            Variable *m13,
            CountingVariable *eventNumber,
            DecayInfo *decay,
            std::vector<MixingTimeResolution *> &r,
            GooPdf *eff,
            Variable *md0,
            Variable *mistag = nullptr);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalisation will get *really* confused and give wrong answers.

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
    // Note that normalisation is not affected because the integrals of S and S'
    // are identical and the weights sum to one, and efficiency is not affected
    // because it depends on the momenta of the daughter tracks, which are not
    // affected by making the wrong charge assignment to the mother.

    __host__ fptype normalize() const override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 5);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }

  protected:
  private:
    DecayInfo *decayInfo;
    Variable *_m12;
    Variable *_m13;
    fptype *dalitzNormRange{nullptr};

    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<WaveHolder_s> *cachedWaves[16]; // Caches the BW values for each event.
    ThreeComplex ***integrals{nullptr}; // Caches the integrals of the BW waves for each combination of resonances.

    bool *redoIntegral;
    mutable bool forceRedoIntegrals{true};
    fptype *cachedMasses;
    fptype *cachedWidths;
    MixingTimeResolution *resolution;

    int totalEventSize;
    int cacheToUse{0};
    SpecialDalitzIntegrator ***integrators{nullptr};
    SpecialWaveCalculator **calculators{nullptr};
};

class SpecialDalitzIntegrator : public thrust::unary_function<thrust::tuple<int, fptype *>, ThreeComplex> {
  public:
    SpecialDalitzIntegrator(int pIdx, unsigned int ri, unsigned int rj);
    __device__ ThreeComplex operator()(thrust::tuple<int, fptype *> t) const;

  private:
    unsigned int resonance_i;
    unsigned int resonance_j;
    unsigned int parameters;
};

class SpecialComplexSum : public thrust::binary_function<ThreeComplex, ThreeComplex, ThreeComplex> {
  public:
    __host__ __device__ ThreeComplex operator()(ThreeComplex one, ThreeComplex two) {
        return ThreeComplex(thrust::get<0>(one) + thrust::get<0>(two),
                            thrust::get<1>(one) + thrust::get<1>(two),
                            thrust::get<2>(one) + thrust::get<2>(two),
                            thrust::get<3>(one) + thrust::get<3>(two),
                            thrust::get<4>(one) + thrust::get<4>(two),
                            thrust::get<5>(one) + thrust::get<5>(two));
    }
};

class SpecialWaveCalculator : public thrust::unary_function<thrust::tuple<int, fptype *, int>, WaveHolder_s> {
  public:
    SpecialWaveCalculator(int pIdx, unsigned int res_idx);
    __device__ WaveHolder_s operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int resonance_i;
    unsigned int parameters;
};

} // namespace GooFit
