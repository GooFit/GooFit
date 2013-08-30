#ifndef FUNCTOR_BASE_HH
#define FUNCTOR_BASE_HH

#include "Variable.hh" 
#include "GlobalCudaDefines.hh" 
#include "FitControl.hh" 
#include <set> 
#include "BinnedDataSet.hh" 
#include "UnbinnedDataSet.hh" 

#include "thrust/iterator/constant_iterator.h" 
#include <thrust/device_vector.h>

typedef thrust::counting_iterator<int> IndexIterator;
typedef thrust::constant_iterator<fptype*> DataIterator;
typedef thrust::constant_iterator<int> SizeIterator;
typedef thrust::tuple<IndexIterator, DataIterator, SizeIterator> EventTuple; 
typedef thrust::zip_iterator<EventTuple> EventIterator; 

const int maxParams = 2000; 
extern fptype* cudaDataArray;
extern fptype host_normalisation[maxParams];
extern fptype host_params[maxParams];
extern unsigned int host_indices[maxParams]; 
extern int totalParams; 
extern int totalConstants;
#ifdef OMP_ON
#pragma omp threadprivate (cudaDataArray, totalParams, totalConstants)
#pragma omp threadprivate (host_normalisation, host_params, host_indices)
#pragma omp threadprivate (host_callnumber)
#endif

class FunctorBase {

public:
  FunctorBase (Variable* x, std::string n); 

  enum Specials {ForceSeparateNorm = 1, ForceCommonNorm = 2}; 

  __host__ virtual double calculateNLL () const = 0; 
  __host__ virtual fptype normalise () const = 0;
  __host__ void initialiseIndices (std::vector<unsigned int> pindices); 

  typedef std::vector<Variable*> obsCont;
  typedef obsCont::iterator obsIter;
  typedef obsCont::const_iterator obsConstIter;
  typedef std::vector<Variable*> parCont; 
  typedef parCont::iterator parIter; 
  typedef parCont::const_iterator parConstIter; 

  __host__ void addSpecialMask (int m) {specialMask |= m;}
  __host__ void copyParams (const std::vector<double>& pars) const;
  __host__ void copyParams (); 
  __host__ void copyNormFactors () const;
  __host__ void generateNormRange (); 
  __host__ std::string getName () const {return name;} 
  __host__ virtual void getObservables (obsCont& ret) const; 
  __host__ virtual void getParameters (parCont& ret) const;
  __host__ Variable* getParameterByName (string n) const; 
  __host__ int getSpecialMask () const {return specialMask;}
  __host__ void setData (BinnedDataSet* data);
  __host__ void setData (UnbinnedDataSet* data);
  __host__ void setData (std::vector<std::map<Variable*, fptype> >& data); 
  __host__ virtual void setFitControl (FitControl* const fc, bool takeOwnerShip = true) = 0; 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 
  __host__ unsigned int getFunctionIndex () const {return functionIdx;}
  __host__ unsigned int getParameterIndex () const {return parameters;}
  __host__ unsigned int registerParameter (Variable* var); 
  __host__ unsigned int registerConstants (unsigned int amount); 
  __host__ virtual void recursiveSetNormalisation (fptype norm = 1) const;
  __host__ void unregisterParameter (Variable* var);
  __host__ void registerObservable (Variable* obs); 
  __host__ void setIntegrationFineness (int i); 
  __host__ void printProfileInfo (bool topLevel = true);

  __host__ bool parametersChanged () const; 
  __host__ void storeParameters () const;

  __host__ obsIter obsBegin () {return observables.begin();}
  __host__ obsIter obsEnd   () {return observables.end();} 
  __host__ obsConstIter obsCBegin () const {return observables.begin();}
  __host__ obsConstIter obsCEnd   () const {return observables.end();} 

  __host__ void checkInitStatus (std::vector<std::string>& unInited) const; 
  void clearCurrentFit (); 

protected:
  fptype numEvents;         // Non-integer to allow weighted events
  unsigned int numEntries;  // Eg number of bins - not always the same as number of events, although it can be. 
  fptype* normRanges;       // This is specific to functor instead of variable so that MetricTaker::operator needn't use indices. 
  unsigned int parameters;  // Stores index, in 'paramIndices', where this functor's information begins. 
  unsigned int cIndex;      // Stores location of constants. 
  obsCont observables; 
  parCont parameterList; 
  FitControl* fitControl; 
  std::vector<FunctorBase*> components;
  int integrationBins; 
  int specialMask; // For storing information unique to PDFs, eg "Normalise me separately" for TddpThrustFunctor. 
  mutable fptype* cachedParams; 
  bool properlyInitialised; // Allows checking for required extra steps in, eg, Tddp and Convolution. 

  unsigned int functionIdx; // Stores index of device function pointer. 

private:
  std::string name; 

  __host__ void recursiveSetIndices ();
  __host__ void setIndices ();
};


#endif
