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
                unsigned int MCeventsNorm = 5e6);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    __host__ fptype normalize() override;

    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 8);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }
    __host__ int getMCevents() { return MCevents; }
    __host__ void setGenerationOffset(int off) { generation_offset = off; }
    __host__ void setMaxWeight(fptype wmax) { maxWeight = wmax; }
    
    __host__ void set_use_bdt_weights(bool use = true){use_bdt_weights = use;}
    
    //functions to get the normalization variables
    //this means copying to host from device and vice-versa but there doesn't seem to be an alternative
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

    __host__ void set_norm_dtime(mcbooster::RealVector_h norm_dtime_h){
      norm_dtime = norm_dtime_h;
    }
       
    __host__ void set_norm_eff(mcbooster::RealVector_h norm_eff_h){
      norm_eff = norm_eff_h;
    }

    __host__ std::
        tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
        GenerateSig(unsigned int numEvents, int seed = 0);

    /*
    __host__ std::
      tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
      GenerateNormEventsWithTime(unsigned int numEvents, int seed = 0);
    */
    __host__ void populateArrays() override;

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
    mcbooster::RealVector_d norm_eff;


    /*
    mcbooster::RealVector_d norm_m12_with_time;
    mcbooster::RealVector_d norm_m34_with_time;
    mcbooster::RealVector_d norm_costheta12_with_time;
    mcbooster::RealVector_d norm_costheta34_with_time;
    mcbooster::RealVector_d norm_phi_with_time;
    */
    
    
    // store spin and lineshape values for normalization
    mutable mcbooster::RealVector_d norm_SF;
    mutable mcbooster::mc_device_vector<fpcomplex> norm_LS;
    mutable bool use_bdt_weights{true};
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
