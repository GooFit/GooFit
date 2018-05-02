/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <mcbooster/GContainers.h>
#include <thrust/remove.h>
#include <tuple>

namespace GooFit {

class LSCalculator_TD;
class AmpCalc_TD;
class SFCalculator_TD;
class NormIntegrator_TD;

class TDDP4 : public GooPdf {
  public:
    TDDP4(std::string n,
          std::vector<Observable> observables,
          DecayInfo4t decay,
          MixingTimeResolution *r,
          GooPdf *eff,
          Observable *mistag,
          unsigned int MCeventsNorm = 5e6);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalisation will get *really* confused and give wrong answers.

    __host__ fptype normalize() const override;

    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 8);

    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }

    __host__ int getMCevents() { return MCevents; }

    __host__ void setGenerationOffset(int off) { generation_offset = off; }

    __host__ void setMaxWeight(fptype wmax) { maxWeight = wmax; }

    __host__ std::
        tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
        GenerateSig(unsigned int numEvents);

    __host__ void recursiveSetIndices() override;

  protected:
  private:
    std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int>>> AmpMap;
    std::map<std::string, unsigned int> compMap;
    // std::map<unsigned int, unsigned int> massmap;
    std::map<std::string, unsigned int> SpinMap;
    std::vector<SpinFactor *> SpinFactors;
    std::vector<Lineshape *> LineShapes;
    std::vector<AmpCalc_TD *> AmpCalcs;
    NormIntegrator_TD *Integrator;
    std::vector<SFCalculator_TD *> sfcalculators;
    std::vector<LSCalculator_TD *> lscalculators;

    unsigned int efficiencyFunction;

    // store normalization events
    mcbooster::RealVector_d norm_M12;
    mcbooster::RealVector_d norm_M34;
    mcbooster::RealVector_d norm_CosTheta12;
    mcbooster::RealVector_d norm_CosTheta34;
    mcbooster::RealVector_d norm_phi;
    // store spin and lineshape values for normalization
    mutable mcbooster::RealVector_d norm_SF;
    mutable mcbooster::mc_device_vector<fpcomplex> norm_LS;

    DecayInfo4t decayInfo;
    MixingTimeResolution *resolution;
    int MCevents;
    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<fpcomplex> *cachedResSF{nullptr}; // Caches the BW values and Spins for each event.
    thrust::device_vector<fpcomplex> *cachedAMPs{nullptr};  // cache Amplitude values for each event.
    mutable bool generation_no_norm{false};
    mutable bool SpinsCalculated{false};
    bool *redoIntegral;
    mutable bool forceRedoIntegrals{true};
    fptype *cachedMasses;
    fptype *cachedWidths;
    int totalEventSize;
    int cacheToUse{0};
    unsigned int generation_offset{25031992};
    double maxWeight{0};
};

class SFCalculator_TD : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fpcomplex> {
  public:
    // Used to create the cached BW values.
    SFCalculator_TD();

    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setSpinFactorId(int idx) { _spinfactor_i = idx; }
    __device__ fpcomplex operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int _spinfactor_i{0};
    unsigned int dalitzFuncId;
};

class NormSpinCalculator_TD
    : public thrust::unary_function<thrust::tuple<fptype, fptype, fptype, fptype, fptype>, fptype> {
  public:
    // Used to create the cached BW values.
    NormSpinCalculator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setSpinFactorId(int idx) { _spinfactor_i = idx; }
    __device__ fptype operator()(thrust::tuple<fptype, fptype, fptype, fptype, fptype> t) const;

  private:
    unsigned int _spinfactor_i{0};
    unsigned int dalitzFuncId;
};

class LSCalculator_TD : public thrust::unary_function<thrust::tuple<int, fptype *, int>, fpcomplex> {
  public:
    // Used to create the cached BW values.
    LSCalculator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setResonanceId(int idx) { _resonance_i = idx; }
    __device__ fpcomplex operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int _resonance_i{0};
    unsigned int dalitzFuncId;
};

class NormLSCalculator_TD
    : public thrust::unary_function<
          thrust::
              tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t>,
          fpcomplex> {
  public:
    // Used to create the cached BW values.
    NormLSCalculator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    void setResonanceId(int idx) { _resonance_i = idx; }
    __device__ fpcomplex operator()(
        thrust::
            tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
        const;

  private:
    unsigned int _resonance_i{0};
    unsigned int dalitzFuncId;
};

class AmpCalc_TD : public thrust::unary_function<unsigned int, fpcomplex> {
  public:
    AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx);
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    // void setAmplitudeId(int idx) { _AmpIdx = idx; }
    // void setpIdx(unsigned int pIdx){_parameters = pIdx;}
    __device__ fpcomplex operator()(thrust::tuple<int, fptype *, int> t) const;

  private:
    unsigned int _nPerm;
    unsigned int _AmpIdx;
    unsigned int dalitzFuncId;
};

class NormIntegrator_TD : public thrust::unary_function<thrust::tuple<int, int, fptype *, fpcomplex *>, fptype> {
  public:
    NormIntegrator_TD();
    void setDalitzId(int idx) { dalitzFuncId = idx; }
    __device__ thrust::tuple<fptype, fptype, fptype, fptype>
    operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const;

  private:
    unsigned int dalitzFuncId;
};

class FourDblTupleAdd : public thrust::binary_function<thrust::tuple<fptype, fptype, fptype, fptype>,
                                                       thrust::tuple<fptype, fptype, fptype, fptype>,
                                                       thrust::tuple<fptype, fptype, fptype, fptype>> {
  public:
    __host__ __device__ thrust::tuple<fptype, fptype, fptype, fptype>
    operator()(thrust::tuple<fptype, fptype, fptype, fptype> one, thrust::tuple<fptype, fptype, fptype, fptype> two) {
        return {thrust::get<0>(one) + thrust::get<0>(two),
                thrust::get<1>(one) + thrust::get<1>(two),
                thrust::get<2>(one) + thrust::get<2>(two),
                thrust::get<3>(one) + thrust::get<3>(two)};
    }
};
} // namespace GooFit
