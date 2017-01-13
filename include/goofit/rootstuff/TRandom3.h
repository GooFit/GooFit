// This file is derived from TRandom3.h in ROOT version 5.30/06.
// It is distributed under the GNU LGPL. The original copyright notice
// reads:

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRandom3
#define ROOT_TRandom3



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRandom3                                                             //
//                                                                      //
// random number generator class: Mersenne Twistor                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TRandom
#include "goofit/rootstuff/TRandom.h"
#endif

class TRandom3 : public TRandom {

private:
    unsigned int   fMt[624];
    int    fCount624;

public:
    TRandom3(unsigned int seed=4357);
    virtual ~TRandom3();
    // get the current seed (only first element of the seed table)
    virtual  unsigned int GetSeed() const {
        return fMt[0];
    }
    virtual  double       Rndm(int i=0);
    virtual  void         RndmArray(int n, float* array);
    virtual  void         RndmArray(int n, double* array);
    virtual  void         SetSeed(unsigned int seed=0);
};

#endif
