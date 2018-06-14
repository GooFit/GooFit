/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

This file includes some lineshapes and spinfactors.
Also right now it is the home to some helper functions needed and an implementation of a simple 4-vec class that works
on the GPU
*/

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

namespace GooFit {

// This constructor is protected
Lineshape::Lineshape(std::string name, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : GooPdf(name)
    , _L(L)
    , _Mpair(Mpair)
    , _FormFac(FormFac)
    , _radius(radius) {
    // Making room for index of decay-related constants. Assumption:
    // These are mother mass and three daughter masses in that order.
    // They will be registered by the object that uses this resonance,
    // which will tell this object where to find them by calling setConstantIndex.
}

} // namespace GooFit
