#pragma once

#include <mcbooster/GContainers.h>

#include <goofit/detail/Complex.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

#include<fstream>

namespace GooFit {

class NormEvents_4Body_DeviceCached final : public NormEvents_4Body_Base {
  public:
    NormEvents_4Body_DeviceCached()                                                       = delete;
    NormEvents_4Body_DeviceCached(const NormEvents_4Body_DeviceCached &copyMe)            = default;
    NormEvents_4Body_DeviceCached(NormEvents_4Body_DeviceCached &&moveMe)                 = default;
    virtual ~NormEvents_4Body_DeviceCached() override                                     = default;
    NormEvents_4Body_DeviceCached &operator=(const NormEvents_4Body_DeviceCached &copyMe) = default;
    NormEvents_4Body_DeviceCached &operator=(NormEvents_4Body_DeviceCached &&moveMe)      = default;

    static std::vector<NormEvents_4Body_Base *>
    buildBatches(const std::vector<long> &normSeeds,
                 unsigned int numNormEventsToGenPerBatch,
                 const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses,
                 bool with_acceptance = false);

    __host__ NormEvents_4Body_DeviceCached(const std::vector<mcbooster::GReal_t> &motherAndDaughterMasses,
                                           long normSeed,
                                           unsigned int numNormEventsToGen);

    __host__ fptype virtual computeNorm_TD(bool noCachedNormValuesToCompute,
                                           const MixingTimeResolution *const resolution,
                                           fptype tau,
                                           fptype xmixing,
                                           fptype ymixing,
                                           unsigned int dalitzId,
                                           bool spinsCalculated,
                                           const std::vector<bool> &lineshapeChanged,
                                           const std::vector<unsigned int> &sfFunctionIndices,
                                           const std::vector<unsigned int> &lsFunctionIndices,
                                           unsigned int CacheIdx) override;
    
    __host__ void copy_norm_phsp_to_host(){
         _norm_M12_h = _norm_M12_d;
         _norm_M34_h = _norm_M34_d;
         _norm_CosTheta12_h = _norm_CosTheta12_d;
         _norm_CosTheta34_h = _norm_CosTheta34_d;
         _norm_phi_h = _norm_phi_d;
      //return std::vector<mcbooster::RealVector_h>{_norm_M12_h,_norm_M34_h,_norm_CosTheta12_h, _norm_CosTheta34_h, _norm_phi_h };
    }

  //Setter function to set normalisation information ,weights, decay time etc.
    __host__ void set_norm_info(mcbooster::RealVector_h norm_dtime_h, mcbooster::RealVector_h norm_dtime_weights_h, mcbooster::RealVector_h norm_eff_h){
      _norm_dtime_d = norm_dtime_h;
      _norm_dtime_weights_d = norm_dtime_weights_h;
      _norm_eff_d = norm_eff_h;
    }

  __host__ std::vector<mcbooster::RealVector_h> get_norm_info(){
    mcbooster::RealVector_h norm_dtime_h = _norm_dtime_d;
    mcbooster::RealVector_h norm_dtime_weights_h = _norm_dtime_weights_d;
    mcbooster::RealVector_h norm_eff_h = _norm_eff_d;
    return std::vector<mcbooster::RealVector_h>{norm_dtime_h,norm_dtime_weights_h, norm_eff_h};
  }

  __host__ int write_norm_phsp(std::string fname){
      std::ofstream outputfile;
      outputfile.open(fname);
      copy_norm_phsp_to_host();
      for(int i = 0; i < _norm_M12_h.size();i++){
        outputfile << _norm_M12_h[i] << " " << _norm_M34_h[i] << " " << _norm_CosTheta12_h[i] << " " << _norm_CosTheta34_h[i] << " " << _norm_phi_h[i] << std::endl;
      }
      outputfile.close();
      printf("Wrote normalisation events");
      return 0;
  }

  __host__ int read_norm_info(std::string fname){
      std::ifstream file(fname);
      std::string line;
      double norm_eff_weight;
      double norm_dtime;
      double norm_dtime_weight;
      mcbooster::RealVector_h norm_dtime_h;
      mcbooster::RealVector_h norm_dtime_weights_h;
      mcbooster::RealVector_h norm_eff_weights_h;
      while (std::getline(file,line)) {
        std::stringstream ss(line);
        if(file >> norm_eff_weight >> norm_dtime >>  norm_dtime_weight){
          norm_dtime_h.push_back(norm_dtime);
          norm_dtime_weights_h.push_back(norm_dtime_weight);
          norm_eff_weights_h.push_back(norm_eff_weight);
        }
      }
      file.close();
      //copy host vectors to device
      _norm_dtime_d = norm_dtime_h;
      _norm_eff_d = norm_eff_weights_h;
      _norm_dtime_weights_d = norm_dtime_weights_h;
    return 0;
  }
  protected:
  private:
    // store normalization events
    mcbooster::RealVector_d _norm_M12_d;
    mcbooster::RealVector_d _norm_M34_d;
    mcbooster::RealVector_d _norm_CosTheta12_d;
    mcbooster::RealVector_d _norm_CosTheta34_d;
    mcbooster::RealVector_d _norm_phi_d;
    // store spin and lineshape values for normalization
    mcbooster::RealVector_d _norm_SF_d;
    mcbooster::mc_device_vector<fpcomplex> _norm_LS_d;
    // efficiency, decay time and importance sample weights
    mcbooster::RealVector_d _norm_dtime_d;
    mcbooster::RealVector_d _norm_eff_d;
    mcbooster::RealVector_d _norm_dtime_weights_d;

    //host vectors of phsp
    mcbooster::RealVector_h _norm_M12_h;
    mcbooster::RealVector_h _norm_M34_h;
    mcbooster::RealVector_h _norm_CosTheta12_h;
    mcbooster::RealVector_h _norm_CosTheta34_h;
    mcbooster::RealVector_h _norm_phi_h;
};

} // end namespace GooFit
