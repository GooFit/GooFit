// This file is derived from TRandom3.cxx in ROOT version 5.30/06.
// It is distributed under the GNU LGPL. The original copyright notice
// is reproduced in full below. 

//////////////////////////////////////////////////////////////////////////
//
// TRandom3
//
// Random number generator class based on
//   M. Matsumoto and T. Nishimura,
//   Mersenne Twistor: A 623-diminsionally equidistributed
//   uniform pseudorandom number generator
//   ACM Transactions on Modeling and Computer Simulation,
//   Vol. 8, No. 1, January 1998, pp 3--30.
//
// For more information see the Mersenne Twistor homepage
//   http://www.math.keio.ac.jp/~matumoto/emt.html
//
// Advantage: large period 2**19937-1
//            relativly fast
//              (only two times slower than TRandom, but
//               two times faster than TRandom2)
// Drawback:  a relative large internal state of 624 integers
//
//
// Aug.99 ROOT implementation based on CLHEP by P.Malzacher
//
// the original code contains the following copyright notice:
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */
/////////////////////////////////////////////////////////////////////

#include "goofit/rootstuff/TRandom3.h"

TRandom* TRandom::gRandom = new TRandom3();

//______________________________________________________________________________
TRandom3::TRandom3 (unsigned int seed) {
//*-*-*-*-*-*-*-*-*-*-*default constructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// If seed is 0, the seed is automatically computed via a TUUID object.
// In this case the seed is guaranteed to be unique in space and time.
   SetSeed(seed);
}

//______________________________________________________________________________
TRandom3::~TRandom3() {
//*-*-*-*-*-*-*-*-*-*-*default destructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ==================
}

//______________________________________________________________________________
double TRandom3::Rndm (int) {
//  Machine independent random number generator.
//  Produces uniformly-distributed floating points in ]0,1]
//  Method: Mersenne Twistor

   unsigned int y;

   const int  kM = 397;
   const int  kN = 624;
   const unsigned int kTemperingMaskB =  0x9d2c5680;
   const unsigned int kTemperingMaskC =  0xefc60000;
   const unsigned int kUpperMask =       0x80000000;
   const unsigned int kLowerMask =       0x7fffffff;
   const unsigned int kMatrixA =         0x9908b0df;

   if (fCount624 >= kN) {
      register int i;

      for (i=0; i < kN-kM; i++) {
         y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
         fMt[i] = fMt[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      for (   ; i < kN-1    ; i++) {
         y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
         fMt[i] = fMt[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      y = (fMt[kN-1] & kUpperMask) | (fMt[0] & kLowerMask);
      fMt[kN-1] = fMt[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      fCount624 = 0;
   }

   y = fMt[fCount624++];
   y ^=  (y >> 11);
   y ^= ((y << 7 ) & kTemperingMaskB );
   y ^= ((y << 15) & kTemperingMaskC );
   y ^=  (y >> 18);

   if (y) return ( (double) y * 2.3283064365386963e-10); // * Power(2,-32)
   return Rndm();
}

//______________________________________________________________________________
void TRandom3::RndmArray(int n, float *array)
{
  // Return an array of n random numbers uniformly distributed in ]0,1]

  for(int i=0; i<n; i++) array[i]=(float)Rndm();
}

//______________________________________________________________________________
void TRandom3::RndmArray(int n, double *array)
{
  // Return an array of n random numbers uniformly distributed in ]0,1]

   int k = 0;

   unsigned int y;

   const int  kM = 397;
   const int  kN = 624;
   const unsigned int kTemperingMaskB =  0x9d2c5680;
   const unsigned int kTemperingMaskC =  0xefc60000;
   const unsigned int kUpperMask =       0x80000000;
   const unsigned int kLowerMask =       0x7fffffff;
   const unsigned int kMatrixA =         0x9908b0df;

   while (k < n) {
      if (fCount624 >= kN) {
         register int i;

         for (i=0; i < kN-kM; i++) {
            y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
            fMt[i] = fMt[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
         }

         for (   ; i < kN-1    ; i++) {
            y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
            fMt[i] = fMt[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
         }

         y = (fMt[kN-1] & kUpperMask) | (fMt[0] & kLowerMask);
         fMt[kN-1] = fMt[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
         fCount624 = 0;
      }

      y = fMt[fCount624++];
      y ^=  (y >> 11);
      y ^= ((y << 7 ) & kTemperingMaskB );
      y ^= ((y << 15) & kTemperingMaskC );
      y ^=  (y >> 18);

      if (y) {
         array[k] = double( y * 2.3283064365386963e-10); // * Power(2,-32)
         k++;
      }
   }
}

//______________________________________________________________________________
void TRandom3::SetSeed(unsigned int seed) {
  //  Set the random generator sequence
  
  TRandom::SetSeed(seed);
  fCount624 = 624;
  int i,j;
  
  fMt[0] = fSeed;
  j = 1;
  // use multipliers from  Knuth's "Art of Computer Programming" Vol. 2, 3rd Ed. p.106
  for(i=j; i<624; i++) {
    fMt[i] = (1812433253 * ( fMt[i-1]  ^ ( fMt[i-1] >> 30)) + i);
  }
}
