#pragma once

#include <goofit/PDFs/physics/AmpComponent.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/detail/Complex.h>

#define NCHANNELS 5

namespace GooFit {

typedef fpcomplex (*resonance_function_ptr)(fptype, fptype, fptype, ParameterContainer &pc);

__device__ auto getERM(fptype mothermass, fptype mbach, fptype mij ) -> fptype;

__device__ auto calcCovFactor( const fptype erm , const int spin)-> fptype;

__device__ auto h(const fptype &m,const fptype &q)->fptype;

__device__ auto h_prime(const fptype &m0,const fptype &q0)->fptype;

__device__ auto d(const fptype &m0,const fptype &q0)->fptype;

__device__ auto f(const fptype &m, const fptype &m0,const fptype &width , const fptype &q, const fptype &q0)->fptype;

__device__ auto calc_q(fptype s12, fptype m1, fptype m2)-> fptype;

__device__ auto DaugDecayMomResFrame(fptype rMassSq, fptype d1m, fptype d2m) -> fptype;

__device__ auto BachMomResFrame(fptype M, fptype rMassSq, fptype mBach) -> fptype;

__device__ auto BachMomParentFrame(fptype rMassSq, fptype M, fptype mBach) -> fptype ;

__device__ auto BlattWeisskopfPrime(fptype z, unsigned int spin)-> fptype;

__device__ auto twoBodyCMmom(double m, fptype m1, fptype m2) -> fptype;

__device__ auto cFromM(fptype motherMass,
                           fptype daug1Mass,
                           fptype daug2Mass,
                           fptype daug3Mass,
                           fptype m13,
                           fptype m23,
                           fptype m12,
                           unsigned int cyclic_index
                           ) -> fptype;

__device__ auto calcLegendrePoly(fptype cosHel, unsigned int spin) -> fptype;

__device__ auto calcZemachSpinFactor(fptype pProd, fptype legPol, unsigned int spin) -> fptype;

__device__ auto twoBodyCMMothermom(fptype rMassSq, fptype dm, fptype d3m) -> fptype;

__device__ auto dampingFactorSquare(const fptype &cmmom, const int &spin, const fptype &mRadius) -> fptype;

__device__ auto dampingFactorSquareNorm(const fptype &cmmom, const int &spin, const fptype &mRadius) -> fptype;

__device__ auto spinFactor(unsigned int spin,
                           fptype motherMass,
                           fptype daug1Mass,
                           fptype daug2Mass,
                           fptype daug3Mass,
                           fptype m12,
                           fptype m13,
                           fptype m23,
                           unsigned int cyclic_index) -> fptype;

__device__ auto phsp_twoBody(fptype s, fptype m0, fptype m1) -> fpcomplex;

__device__ auto phsp_fourPi(fptype s) -> fpcomplex;

__device__ void getCofactor(fptype A[NCHANNELS][NCHANNELS], fptype temp[NCHANNELS][NCHANNELS], int p, int q, int n);

__device__ auto determinant(fptype A[NCHANNELS][NCHANNELS], int n) -> fptype;

__device__ void adjoint(fptype A[NCHANNELS][NCHANNELS], fptype adj[NCHANNELS][NCHANNELS]);

__device__ auto inverse(fptype A[NCHANNELS][NCHANNELS], fptype inverse[NCHANNELS][NCHANNELS]) -> bool;

__device__ void getPropagator(const fptype kMatrix[NCHANNELS][NCHANNELS],
                              const fpcomplex phaseSpace[NCHANNELS],
                              fpcomplex F[NCHANNELS][NCHANNELS],
                              fptype adlerTerm);

__device__ fpcomplex _Vc (fptype c);
__device__ fpcomplex sigma (fpcomplex s, fptype m);
__device__ fpcomplex V(fpcomplex x);
__device__ fpcomplex wp(fpcomplex x);
__device__ fptype q(fptype s);
__device__ fpcomplex tlow(fpcomplex s);
__device__ fpcomplex pn( fpcomplex x, fptype n);
__device__ fpcomplex fu(fpcomplex x);
__device__ fpcomplex Jp(fpcomplex s,fptype m);
__device__ fpcomplex tf0(fpcomplex s);
__device__ fpcomplex t00(fptype s);
__device__ fpcomplex S00(fptype s);
__device__ fptype wp2(fptype s);
__device__ fptype Phim(fptype s);
__device__ fptype derivaPhi(fptype s);
__device__ fptype Phi(fptype s);
__device__ fptype argument(fptype s);
__device__ fptype Inela(fptype s);
__device__ fpcomplex ampt00 (fptype s);


}