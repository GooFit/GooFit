#ifndef TDDP_THRUST_FUNCTOR_HH
#define TDDP_THRUST_FUNCTOR_HH

#include "ThrustPdfFunctor.hh" 
#include "MixingTimeResolution_Aux.hh" 

enum ResaonceType {FLAVOUR = 0, CP_EVEN, CP_ODD};
enum EvaluateType {
  RBW = 0, 
  NONRES, 
  GOU_SAK, 
  GAUSSIAN // Notice that Gaussian uses mass, not mass-squared!
};
enum DaughterPair {PAIR_12 = 0, PAIR_13, PAIR_23}; 

struct ResonanceInfo {
  Variable* amp_real;
  Variable* amp_imag;
  Variable* mass;
  Variable* width;
  unsigned int spin;
  unsigned int cyclic_index;
  unsigned int eval_type;
  unsigned int resonance_type; 
};

struct DecayInfo {
  fptype motherMass;
  fptype daug1Mass;
  fptype daug2Mass;
  fptype daug3Mass;
  fptype meson_radius;

  Variable* _tau; 
  Variable* _xmixing;
  Variable* _ymixing;
  std::vector<ResonanceInfo*> resonances; 
};

typedef thrust::tuple<fptype, fptype, fptype, fptype> WaveHolder; 
typedef thrust::tuple<fptype, fptype, fptype, fptype, fptype, fptype> ThreeComplex;

class SpecialDalitzIntegrator;
class SpecialWaveCalculator; 

// Copied from strided_range thrust example by Nathan Bell.
// Iterator to move forward by a specified number of steps 
// in each iteration.
template <typename Iterator> class strided_range { 
public:
  typedef typename thrust::iterator_difference<Iterator>::type difference_type;

  struct stride_functor : public thrust::unary_function<difference_type,difference_type> {
    difference_type stride;
    
    stride_functor(difference_type stride)
      : stride(stride) {}
    
    __host__ __device__ difference_type operator() (const difference_type& i) const {
      return stride * i;
    }
  };
  typedef typename thrust::counting_iterator<difference_type>                   CountingIterator;
  typedef typename thrust::transform_iterator<stride_functor, CountingIterator> TransformIterator;
  typedef typename thrust::permutation_iterator<Iterator,TransformIterator>     PermutationIterator;

  // type of the strided_range iterator
  typedef PermutationIterator iterator;

  // construct strided_range for the range [first,last)
  strided_range (Iterator first, Iterator last, difference_type stride)
    : first(first), last(last), stride(stride) {}
   
  iterator begin () const {
    return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride)));
  }

  iterator end () const {
    return begin() + ((last - first) + (stride - 1)) / stride;
  }
  
protected:
  Iterator first;
  Iterator last;
  difference_type stride;
};
  
class TddpThrustFunctor : public ThrustPdfFunctor {
public:
  TddpThrustFunctor (std::string n, Variable* _dtime, Variable* _sigmat, Variable* m12, Variable* m13, Variable* eventNumber, DecayInfo* decay, MixingTimeResolution* r, ThrustPdfFunctor* eff);
  // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the 
  // coherent sum. The caching method requires that it be done this way or the ProdThrustFunctor
  // normalisation will get *really* confused and give wrong answers. 

  __host__ virtual fptype normalise () const;
  __host__ void setDataSize (unsigned int dataSize, unsigned int evtSize = 5); 
  __host__ void setForceIntegrals (bool f = true) {forceRedoIntegrals = f;}  

protected:

private:
  DecayInfo* decayInfo; 
  Variable* _m12;
  Variable* _m13; 
  fptype* dalitzNormRange; 

  // Following variables are useful if masses and widths, involved in difficult BW calculation, 
  // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
  thrust::device_vector<WaveHolder>* cachedWaves; // Caches the BW values for each event.
  ThreeComplex*** integrals; // Caches the integrals of the BW waves for each combination of resonances. 

  bool* redoIntegral;
  mutable bool forceRedoIntegrals; 
  fptype* cachedMasses; 
  fptype* cachedWidths;
  MixingTimeResolution* resolution; 
  int totalEventSize; 
  int cacheToUse; 
  SpecialDalitzIntegrator*** integrators;
  SpecialWaveCalculator** calculators; 
};

class SpecialDalitzIntegrator : public thrust::unary_function<thrust::tuple<int, fptype*>, ThreeComplex > {
public:

  SpecialDalitzIntegrator (int pIdx, unsigned int ri, unsigned int rj);
  __device__ ThreeComplex operator () (thrust::tuple<int, fptype*> t) const;
private:

  unsigned int resonance_i;
  unsigned int resonance_j; 
  unsigned int parameters;
}; 

class SpecialComplexSum : public thrust::binary_function<ThreeComplex,
							 ThreeComplex,
							 ThreeComplex > {

public: 

  __host__ __device__ ThreeComplex operator() (ThreeComplex one, ThreeComplex two) {
    return ThreeComplex(thrust::get<0>(one) + thrust::get<0>(two), 
			thrust::get<1>(one) + thrust::get<1>(two), 
			thrust::get<2>(one) + thrust::get<2>(two), 
			thrust::get<3>(one) + thrust::get<3>(two), 
			thrust::get<4>(one) + thrust::get<4>(two), 
			thrust::get<5>(one) + thrust::get<5>(two));
  }
};

class SpecialWaveCalculator : public thrust::unary_function<thrust::tuple<int, fptype*, int>, WaveHolder > {
public:

  SpecialWaveCalculator (int pIdx, unsigned int res_idx); 
  __device__ WaveHolder operator () (thrust::tuple<int, fptype*, int> t) const;

private:

  unsigned int resonance_i;
  unsigned int parameters;
}; 


#endif

