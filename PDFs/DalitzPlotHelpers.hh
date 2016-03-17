#ifndef DALITZ_PLOT_HELPERS_HH
#define DALITZ_PLOT_HELPERS_HH

#include "ResonancePdf.hh"
#include "LineshapesPdf.hh"



EXEC_TARGET bool inDalitz (fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3); 
EXEC_TARGET devcomplex<fptype> getResonanceAmplitude (fptype m12, fptype m13, fptype m23, 
						      unsigned int functionIdx, unsigned int pIndex); 

EXEC_TARGET void get4Vecs (fptype* Vecs, const unsigned int& constants, const fptype& m12, const fptype& m34, const fptype& cos12, const fptype& cos34, const fptype& phi){
  fptype M = functorConstants[constants + 1]; 
  fptype m1  = functorConstants[constants + 2]; 
  fptype m2  = functorConstants[constants + 3]; 
  fptype m3  = functorConstants[constants + 4]; 
  fptype m4  = functorConstants[constants + 5]; 
  
  fptype E1 = (m12*m12 + m1*m1 - m2*m2) / (2 * m12) ; 
  fptype E2 = (m12*m12 - m1*m1 + m2*m2) / (2 * m12) ; 
  fptype E3 = (m34*m34 + m3*m3 - m4*m4) / (2 * m34) ; 
  fptype E4 = (m34*m34 - m3*m3 + m4*m4) / (2 * m34) ; 
  fptype p1 = sqrt(E1*E1 - m1*m1);
  fptype p3 = sqrt(E3*E3 - m3*m3); 
  fptype sin12 = sqrt(1-cos12*cos12);
  fptype sin34 = sqrt(1-cos34*cos34);
  fptype ED1 = ( M*M + m12*m12 - m34*m34) / (2*m12);
  fptype PD1 = sqrt(ED1*ED1 - M*M);
  fptype beta1 = PD1 / ED1;
  fptype gamma1 = 1.0/sqrt(1-beta1*beta1);
  fptype ED2 = ( M*M - m12*m12 + m34*m34) / (2*m34);
  fptype PD2 = sqrt(ED2*ED2 - M*M);
  fptype beta2 = -PD2 / ED2;
  fptype gamma2 = 1.0/sqrt(1-beta2*beta2);

  //set X-component
  Vecs[0] = cos12*p1;
  Vecs[4] = -cos12*p1;
  Vecs[8] = -cos34*p3;
  Vecs[12] = cos34*p3;
 
  //set Y-component
  Vecs[1] = sin12*p1;
  Vecs[5] = -sin12*p1;
  Vecs[9] = -sin34*p3;
  Vecs[13] = sin34*p3;

  //set E-component
  Vecs[3]  = E1;
  Vecs[7]  = E2;
  Vecs[11] = E3;
  Vecs[15] = E4;

  fptype tmpE = Vecs[3];
  fptype tmpX = Vecs[0];
  Vecs[3] = gamma1*( tmpE + beta1*tmpX );
  Vecs[0] = gamma1*( tmpX + beta1*tmpE );

  tmpE = Vecs[7];
   tmpX = Vecs[4];
  Vecs[7] = gamma1*( tmpE + beta1*tmpX );
  Vecs[4] = gamma1*( tmpX + beta1*tmpE );

  tmpE = Vecs[11];
  tmpX = Vecs[8];
  Vecs[11] = gamma2*( tmpE + beta2*tmpX );
  Vecs[8] = gamma2*( tmpX + beta2*tmpE );

  tmpE = Vecs[15];
  tmpX = Vecs[12];
  Vecs[15] = gamma2*( tmpE + beta2*tmpX );
  Vecs[12] = gamma2*( tmpX + beta2*tmpE );

  // rotation around X-axis of the first two vectors.
  fptype cosphi = cos(phi);
  fptype sinphi = sin(phi);

  // note that Z-component is zero thus rotation is as easy as this:
  Vecs[2] = sinphi*Vecs[1]; 
  Vecs[1] = cosphi*Vecs[1];

  Vecs[6] = sinphi*Vecs[5]; 
  Vecs[5] = cosphi*Vecs[5];
  
} 


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
  std::map<unsigned int, fptype> particle_masses ;
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
