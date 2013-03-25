// This file is derived from TRandom.h in ROOT version 5.30/06.
// It is distributed under the GNU LGPL. The original copyright notice
// reads:

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRandom
#define ROOT_TRandom

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRandom                                                              //
//                                                                      //
// Simple prototype random number generator class (periodicity = 10**9) //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TRandom {

protected:
   unsigned int   fSeed;  //Random number generator seed

public:
   TRandom(unsigned int seed=65539);
   virtual ~TRandom();
   virtual  int    Binomial(int ntot, double prob);
   virtual  double BreitWigner(double mean=0, double gamma=1);
   virtual  void     Circle(double &x, double &y, double r);
   virtual  double Exp(double tau);
   virtual  double Gaus(double mean=0, double sigma=1);
   virtual  unsigned int   GetSeed() const {return fSeed;}
   virtual  unsigned int   Integer(unsigned int imax);
   //virtual  double Landau(double mean=0, double sigma=1); // Dropped so I don't have to rip landau_quantile. 
   virtual  int    Poisson(double mean);
   virtual  double PoissonD(double mean);
   virtual  void     Rannor(float &a, float &b);
   virtual  void     Rannor(double &a, double &b);
   virtual  void     SetSeed(unsigned int seed=0);
   virtual  double Rndm(int i=0);
   virtual  void     RndmArray(int n, float *array);
   virtual  void     RndmArray(int n, double *array);
   virtual  void     Sphere(double &x, double &y, double &z, double r);
   virtual  double Uniform(double x1=1);
   virtual  double Uniform(double x1, double x2);

   static TRandom *gRandom;
};

#endif
