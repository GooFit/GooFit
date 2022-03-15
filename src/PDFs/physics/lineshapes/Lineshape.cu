/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficiently tested yet and still under heavy development!

This file includes some lineshapes and spinfactors.
Also right now it is the home to some helper functions needed and an implementation of a simple 4-vec class that works
on the GPU
*/

#include <vector>
#include <string>

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

namespace GooFit {

// This constructor is protected
Lineshape::Lineshape(
    std::string pdf_name, std::string name, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius)
    : AmpComponent("Lineshapes::" + pdf_name, name)
    , _L(L)
    , _Mpair(Mpair)
    , _FormFac(FormFac)
    , _radius(radius) {
    // Making room for index of decay-related constants. Assumption:
    // These are mother mass and three daughter masses in that order.
    // They will be registered by the object that uses this resonance,
    // which will tell this object where to find them by calling setConstantIndex.

    registerConstant(Mpair);
}

bool Lineshape::Lineshape::isEqualByValue(const Lineshape &other) const {
    return this->getName() == other.getName() && this->_L == other._L && this->_Mpair == other._Mpair
           && this->_FormFac == other._FormFac && this->_radius == other._radius
           && this->areParamsandConstantsEqualByVal(other);
}

std::ostream &operator<<(std::ostream &out, const Lineshape &obj) { return obj.print(out); }

std::ostream &Lineshape::Lineshape::print(std::ostream &out) const {
    std::string paramNames = "";
    for(int p = 0; p < this->getParameters().size(); p++) {
        paramNames += this->getParameters()[p].getName() + ", ";
    }

    DP4Pair mpairEnum = static_cast<DP4Pair>(this->_Mpair);

    out << this->getPdfName() << ", " << this->getName() << ", " << paramNames << this->_L << ", " << mpairEnum << ", "
        << this->_FormFac;
    return out;
}

std::ostream &operator<<(std::ostream &out, const FF &obj) {
    std::string ffName;

    switch(obj) {
    case FF::One:
        ffName = "FF::One";
        break;
    case FF::BL:
        ffName = "FF::BL";
        break;
    case FF::BL_Prime:
        ffName = "FF::BL_Prime";
        break;
    case FF::BL2:
        ffName = "FF::BL2";
        break;
    default:
        ffName = "UNKNOWN";
        break;
    }

    out << ffName;
    return out;
}

} // namespace GooFit
