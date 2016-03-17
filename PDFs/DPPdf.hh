#ifndef D_P_PDF_HH
#define D_P_PDF_HH

#include "GooPdf.hh" 
#include "DalitzPlotHelpers.hh" 
#include "devcomplex.hh"

class SpecialIntegrator;
class LSCalculator; 
class AmpCalc;

class DPPdf : public GooPdf {
public:
  DPPdf (std::string n, std::vector<Variable*> observables, Variable* eventNumber, DecayInfo_DP* decay, GooPdf* eff);
  // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the 
  // coherent sum. The caching method requires that it be done this way or the ProdPdf
  // normalisation will get *really* confused and give wrong answers. 

  __host__ virtual fptype normalise () const;
  __host__ void setDataSize (unsigned int dataSize, unsigned int evtSize = 3); 
  __host__ void setForceIntegrals (bool f = true) {forceRedoIntegrals = f;}  
protected:

private:

  std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int> > > AmpMap;
  std::map<std::string, unsigned int> compMap;
  std::map<unsigned int, unsigned int> massmap;
  std::map<std::string, unsigned int> SpinMap;
  std::vector<SpinFactor*> SpinFactors;
  std::vector<AmpCalc*> AmpCalcs;


  DecayInfo_DP* decayInfo; 
  std::vector<Variable*> _observables; 
  fptype* dalitzNormRange; 

  // Following variables are useful if masses and widths, involved in difficult BW calculation, 
  // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
  DEVICE_VECTOR<devcomplex<fptype> >* cachedResSF; // Caches the BW values and Spins for each event.
  DEVICE_VECTOR<devcomplex<fptype> >* cachedAMPs; // cache Amplitude values for each event.
  devcomplex<fptype>*** integrals; // Caches the integrals of the BW waves for each combination of resonances. 

  mutable bool SpinsCalculated;
  bool* redoIntegral;
  mutable bool forceRedoIntegrals; 
  fptype* cachedMasses; 
  fptype* cachedWidths;
  int totalEventSize; 
  int cacheToUse; 
  LSCalculator** calculators;
  SpecialIntegrator*** integrators;
};


class SFCalculator : public thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype> > {
public:
  // Used to create the cached BW values. 
  SFCalculator (int pIdx, unsigned int sf_idx); 
  EXEC_TARGET devcomplex<fptype> operator () (thrust::tuple<int, fptype*, int> t) const;

private:

  unsigned int _spinfactor_i;
  unsigned int _parameters;
}; 


class LSCalculator : public thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype> > {
public:
  // Used to create the cached BW values. 
  LSCalculator (int pIdx, unsigned int res_idx); 
  EXEC_TARGET devcomplex<fptype> operator () (thrust::tuple<int, fptype*, int> t) const;

private:

  unsigned int _resonance_i;
  unsigned int _parameters;
}; 


class AmpCalc : public thrust::unary_function<unsigned int, devcomplex<fptype> >{
  public:
    AmpCalc(unsigned int AmpIdx, unsigned int pIdx, unsigned int CoeffIdx);
    // void setpIdx(unsigned int pIdx){_parameters = pIdx;}
    EXEC_TARGET devcomplex<fptype> operator() (thrust::tuple<int, fptype*, int> t) const;
  private:
    unsigned int _coeffIdx;
    unsigned int _AmpIdx;
    unsigned int _parameters;
 };


class SpecialIntegrator : public thrust::unary_function<thrust::tuple<int, fptype*>, devcomplex<fptype> > {
public:
  // Class used to calculate integrals of terms BW_i * BW_j^*. 
  SpecialIntegrator (int pIdx, unsigned int ri, unsigned int rj, unsigned int starti, unsigned int startj);
  EXEC_TARGET devcomplex<fptype> operator () (thrust::tuple<int, fptype*> t) const;
private:

  unsigned int _amp_i;
  unsigned int _starti;
  unsigned int _startj;
  unsigned int _amp_j; 
  unsigned int _parameters;
}; 

#endif

