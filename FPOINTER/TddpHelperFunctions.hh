#ifndef TDDP_HELPER_HH
#define TDDP_HELPER_HH

#include "devcomplex.hh" 

__device__ fptype twoBodyCMmom (double rMassSq, fptype d1m, fptype d2m);
__device__ fptype dampingFactorSquare (fptype cmmom, int spin, fptype mRadius);
__device__ bool inDalitz (fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3);
__device__ fptype spinFactor (unsigned int spin, fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, fptype m12, fptype m13, fptype m23, unsigned int cyclic_index);
__device__ devcomplex<fptype> plainBW (fptype m12, fptype m13, fptype m23, 
				       fptype resmass, fptype reswidth, fptype meson_radius, 
				       fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, 
				       unsigned int spin, unsigned int cyclic_index); 
__device__ devcomplex<fptype> gaussian (fptype m12, fptype m13, fptype m23, fptype resmass, fptype reswidth, unsigned int cyclic_index);
__device__ void getAmplitudeCoefficients (devcomplex<fptype> a1, devcomplex<fptype> a2, fptype& a1sq, fptype& a2sq, fptype& a1a2real, fptype& a1a2imag);
__device__ devcomplex<fptype> getResonanceAmplitude (fptype m12, fptype m13, fptype resmass, fptype reswidth, unsigned int spin, unsigned int cyclic_index, unsigned int eval_type, unsigned int* indices); 

#endif
