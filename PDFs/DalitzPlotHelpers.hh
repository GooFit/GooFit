/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#ifndef DALITZ_PLOT_HELPERS_HH
#define DALITZ_PLOT_HELPERS_HH

#include "GooPdf.hh" 
#include "ResonancePdf.hh"
#include "LineshapesPdf.hh"

EXEC_TARGET bool inDalitz (fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3); 
EXEC_TARGET devcomplex<fptype> getResonanceAmplitude (fptype m12, fptype m13, fptype m23, unsigned int functionIdx, unsigned int pIndex); 
EXEC_TARGET void get4Vecs (fptype* Vecs, const unsigned int& constants, const fptype& m12, const fptype& m34, const fptype& cos12, const fptype& cos34, const fptype& phi);
EXEC_TARGET fptype getmass(const unsigned int& pair, fptype& d1, fptype& d2, const fptype* vecs, const fptype& m1, const fptype& m2, const fptype& m3, const fptype& m4);

class ALIGN(16) gpuLVec
{
  private:
    fptype X;
    fptype Y;
    fptype Z;
    fptype E;
  public:
    EXEC_TARGET gpuLVec(fptype x, fptype y, fptype z, fptype e);
    EXEC_TARGET fptype getX() const {return X;}
    EXEC_TARGET fptype getY() const {return Y;}
    EXEC_TARGET fptype getZ() const {return Z;}
    EXEC_TARGET fptype getE() const {return E;}
    EXEC_TARGET fptype Dot(const gpuLVec& rhs) const;
    EXEC_TARGET gpuLVec& operator+=(const gpuLVec& rhs);
    EXEC_TARGET gpuLVec& operator-=(const gpuLVec& rhs);
    EXEC_TARGET gpuLVec& operator*=(const fptype rhs);
    // EXEC_TARGET gpuLVec& operator/=(const fptype rhs);
};

EXEC_TARGET gpuLVec operator+(gpuLVec lhs, const gpuLVec& rhs);
EXEC_TARGET gpuLVec operator-(gpuLVec lhs, const gpuLVec& rhs);
EXEC_TARGET gpuLVec operator*(gpuLVec lhs, fptype rhs);
EXEC_TARGET gpuLVec operator*(fptype lhs, gpuLVec rhs);

// in case of 3 particles the first two are the resonance.
enum DP4Pair {M_12 = 0, M_34, M_13, M_14, M_23, M_24, M_123, M_132, M_231, M_124, M_142, M_241, M_134, M_143, M_341, M_234, M_243, M_342}; 
enum DaughterPair {PAIR_12 = 0, PAIR_13, PAIR_23}; 

const int resonanceSize = 4;   // Number of parameters to describe one resonance.
// Why not make this a static const member of ResonancePdf? Because the 'static'
// keyword (and 'extern' as well) interacts badly with some nvcc versions when the
// variable is used in device code. 

struct DecayInfo {
  fptype motherMass;
  fptype daug1Mass;
  fptype daug2Mass;
  fptype daug3Mass;
  fptype meson_radius;

  Variable* _tau; 
  Variable* _xmixing;
  Variable* _ymixing;
  std::vector<ResonancePdf*> resonances; 
};

struct DecayInfo_DP {
  std::vector<fptype> particle_masses ;
  fptype meson_radius;
  std::vector<Amplitude*> amplitudes; 
};

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


#endif
