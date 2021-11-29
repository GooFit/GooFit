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

class Amp4Body_TD : public Amp4BodyBase {
  public:
    Amp4Body_TD(std::string n,
                std::vector<Observable> observables,
                DecayInfo4t decay,
                MixingTimeResolution *r,
                GooPdf *eff,
                Observable *mistag,
                unsigned int MCeventsNorm = 5e6,
                bool special_integration=false,
                unsigned int generationOffset=0);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    __host__ auto normalize() -> fptype override;

    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 8);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }
    __host__ auto getMCevents() -> int { return MCevents; }
    __host__ void setGenerationOffset(int off) { generation_offset = off; }
    __host__ void setMaxWeight(fptype wmax) { maxWeight = wmax; }

    __host__ auto GenerateSig(unsigned int numEvents, int seed = 0) -> std::
        tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>;

    __host__ unsigned int getGenerationOffset(){return generation_offset;}
    __host__ void populateArrays() override;
        __host__ mcbooster::RealVector_h get_norm_m12(){
      auto host_norm_m12 = mcbooster::RealVector_h(norm_M12);
      return host_norm_m12;
    }
     __host__ mcbooster::RealVector_h get_norm_m34()
    {
      auto host_norm_m34 = mcbooster::RealVector_h(norm_M34);
      return host_norm_m34;
    }
    __host__ mcbooster::RealVector_h get_norm_c12(){
      auto host_norm_c12 = mcbooster::RealVector_h(norm_CosTheta12);
      return host_norm_c12;
    }
    __host__ mcbooster::RealVector_h get_norm_c34(){
      auto host_norm_c34 = mcbooster::RealVector_h(norm_CosTheta34);
      return host_norm_c34;
    }
    __host__ mcbooster::RealVector_h get_norm_phi(){
      auto host_norm_phi = mcbooster::RealVector_h(norm_phi);
      return host_norm_phi;
    }

    __host__ mcbooster::RealVector_h get_norm_dtime(){
      auto host_norm_dtime = mcbooster::RealVector_h(norm_dtime);
      return host_norm_dtime;
    }

    __host__ mcbooster::RealVector_h get_norm_eff(){
      auto host_norm_eff = mcbooster::RealVector_h(norm_eff);
      return host_norm_eff;
    }

    __host__ mcbooster::RealVector_h get_norm_importance_weights(){
      auto host_importance_weight = mcbooster::RealVector_h(norm_importance_weight);
      return host_importance_weight;
    }

    __host__ void set_norm_dtime(mcbooster::RealVector_h norm_dtime_h){
      norm_dtime = norm_dtime_h;
    }

    __host__ void set_norm_eff(mcbooster::RealVector_h norm_eff_h){
      norm_eff = norm_eff_h;
    }

    __host__ void set_norm_importance_weights(mcbooster::RealVector_h norm_importance_h){
      norm_importance_weight = norm_importance_h;
    }

     __host__ void set_special_integral(bool special){
      specialIntegral = special;
    }

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
    mcbooster::RealVector_d norm_dtime;
    //efficiency weight given by BDT
    mcbooster::RealVector_d norm_eff;
    //weights from Importance Sampling stays constant throughout 
    mcbooster::RealVector_d norm_importance_weight;
    //determine whether to use full numerical integration for per-event efficiencies
    mutable bool specialIntegral{false};
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
    unsigned int generation_offset{0};
    double maxWeight{0};
};

} // namespace GooFit
