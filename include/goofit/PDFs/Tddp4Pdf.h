/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#ifndef TDDP4_PDF_HH
#define TDDP4_PDF_HH

#include "goofit/PDFs/GooPdf.h" 
#include "goofit/PDFs/DalitzPlotHelpers.h" 
#include "goofit/PDFs/SpinFactors.h"
#include <mcbooster/GContainers.h>
#include <tuple>
#include <thrust/remove.h>
#include "goofit/PDFs/MixingTimeResolution_Aux.h" 
class LSCalculator_TD; 
class AmpCalc_TD;
class SFCalculator_TD;
class NormIntegrator_TD;

class TDDP4 : public GooPdf {
public:
  TDDP4 (std::string n, std::vector<Variable*> observables, DecayInfo_DP* decay, MixingTimeResolution* r, GooPdf* eff, Variable* mistag = 0, unsigned int MCeventsNorm = 5e6);
  // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the 
  // coherent sum. The caching method requires that it be done this way or the ProdPdf
  // normalisation will get *really* confused and give wrong answers. 

  __host__ virtual fptype normalise () const;
  __host__ void setDataSize (unsigned int dataSize, unsigned int evtSize = 8); 
  __host__ void setForceIntegrals (bool f = true) {forceRedoIntegrals = f;}  
  __host__ int getMCevents(){return MCevents;}
  __host__ void setGenerationOffset(int off){generation_offset = off;}
  // __host__ void setGenDecayTimeLimit(double low, double high){genlow = exp(-high); genhigh = exp(-low);}
  __host__ void setGenDecayTimeLimit(double low, double high){genlow = low; genhigh = high;}
  __host__ std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h,  mcbooster::RealVector_h> GenerateSig (unsigned int numEvents);

protected:

private:

  std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int> > > AmpMap;
  std::map<std::string, unsigned int> compMap;
  // std::map<unsigned int, unsigned int> massmap;
  std::map<std::string, unsigned int> SpinMap;
  std::vector<SpinFactor*> SpinFactors;
  std::vector<AmpCalc_TD*> AmpCalcs;
  NormIntegrator_TD* Integrator;
  std::vector<SFCalculator_TD*> sfcalculators;
  std::vector<LSCalculator_TD*> lscalculators;

  // store normalization events
  mcbooster::RealVector_d norm_M12;
  mcbooster::RealVector_d norm_M34;
  mcbooster::RealVector_d norm_CosTheta12;
  mcbooster::RealVector_d norm_CosTheta34;
  mcbooster::RealVector_d norm_phi;
  //store spin and lineshape values for normalization
  mutable mcbooster::RealVector_d norm_SF;
  mutable mcbooster::mc_device_vector<devcomplex<fptype> > norm_LS;

  DecayInfo_DP* decayInfo; 
  std::vector<Variable*> _observables; 
  MixingTimeResolution* resolution; 
  fptype* hostphsp;
  int MCevents;
  // Following variables are useful if masses and widths, involved in difficult BW calculation, 
  // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
  DEVICE_VECTOR<devcomplex<fptype> >* cachedResSF; // Caches the BW values and Spins for each event.
  DEVICE_VECTOR<devcomplex<fptype> >* cachedAMPs; // cache Amplitude values for each event.
  mutable bool generation_no_norm;
  mutable bool SpinsCalculated;
  bool* redoIntegral;
  mutable bool forceRedoIntegrals; 
  fptype* cachedMasses; 
  fptype* cachedWidths;
  int totalEventSize; 
  int cacheToUse; 
  unsigned int generation_offset;
  double genlow;
  double genhigh;
};

class SFCalculator_TD : public thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype> > {
public:
  // Used to create the cached BW values. 
  SFCalculator_TD (int pIdx, unsigned int sf_idx); 
  EXEC_TARGET devcomplex<fptype> operator () (thrust::tuple<int, fptype*, int> t) const;

private:

  unsigned int _spinfactor_i;
  unsigned int _parameters;
}; 

class NormSpinCalculator_TD : public thrust::unary_function<thrust::tuple<fptype, fptype, fptype, fptype, fptype>, fptype> {
public:
  // Used to create the cached BW values. 
  NormSpinCalculator_TD (int pIdx, unsigned int sf_idx); 
  EXEC_TARGET fptype operator () (thrust::tuple<fptype, fptype, fptype, fptype, fptype> t) const;

private:

  unsigned int _spinfactor_i;
  unsigned int _parameters;
}; 


class LSCalculator_TD : public thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype> > {
public:
  // Used to create the cached BW values. 
  LSCalculator_TD (int pIdx, unsigned int res_idx); 
  EXEC_TARGET devcomplex<fptype> operator () (thrust::tuple<int, fptype*, int> t) const;

private:

  unsigned int _resonance_i;
  unsigned int _parameters;
}; 

class NormLSCalculator_TD : public thrust::unary_function<thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t>, devcomplex<fptype> > {
public:
  // Used to create the cached BW values. 
  NormLSCalculator_TD (int pIdx, unsigned int res_idx); 
  EXEC_TARGET devcomplex<fptype> operator () (thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t) const;

private:

  unsigned int _resonance_i;
  unsigned int _parameters;
}; 

class AmpCalc_TD : public thrust::unary_function<unsigned int, devcomplex<fptype> >{
  public:
    AmpCalc_TD(unsigned int AmpIdx, unsigned int pIdx, unsigned int nPerm);
    // void setpIdx(unsigned int pIdx){_parameters = pIdx;}
    EXEC_TARGET devcomplex<fptype> operator() (thrust::tuple<int, fptype*, int> t) const;
  private:
    unsigned int _nPerm;
    unsigned int _AmpIdx;
    unsigned int _parameters;
 };

class NormIntegrator_TD : public thrust::unary_function<thrust::tuple<int, int, fptype*, devcomplex<fptype>*>, fptype >{
  public:
    NormIntegrator_TD(unsigned int pIdx);
    EXEC_TARGET thrust::tuple<fptype,fptype,fptype,fptype> operator() (thrust::tuple<int, int, fptype*, devcomplex<fptype>*> t) const;
  private:
    unsigned int _parameters;
 };

 class FourDblTupleAdd : public thrust::binary_function<thrust::tuple<fptype,fptype,fptype,fptype>,
               thrust::tuple<fptype,fptype,fptype,fptype>,
               thrust::tuple<fptype,fptype,fptype,fptype> > {

 public: 

   __host__ __device__ thrust::tuple<fptype,fptype,fptype,fptype> operator() (thrust::tuple<fptype,fptype,fptype,fptype> one, thrust::tuple<fptype,fptype,fptype,fptype> two) {
     return thrust::tuple<fptype,fptype,fptype,fptype>(thrust::get<0>(one) + thrust::get<0>(two), 
      thrust::get<1>(one) + thrust::get<1>(two),
      thrust::get<2>(one) + thrust::get<2>(two),
      thrust::get<3>(one) + thrust::get<3>(two));
   }
 };

#endif

