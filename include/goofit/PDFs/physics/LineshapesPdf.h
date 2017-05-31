/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/physics/ResonancePdf.h"
#include "goofit/Variable.h"

#include <thrust/complex.h>

namespace GooFit {

class SpinFactor;

enum class LS { ONE, BW, Lass, Lass_M3, nonRes, Bugg, Bugg3, Flatte, SBW };
// PDG notation for FF
enum class FF : unsigned int { One = 0, BL, BL_Prime, BL2 };

class Lineshape : public GooPdf {
    friend class DPPdf;
    friend class TDDP4;
    // Service class intended to hold parametrisations of
    // resonances on Dalitz plots. Don't try to use this
    // as a standalone PDF! It should only be used as a
    // component in one of the friend classes. It extends
    // GooPdf so as to take advantage of the
    // infrastructure, but will crash if used on its own.
    Variable *_mass;
    Variable *_width;
    unsigned int _L;
    unsigned int _Mpair;
    LS _kind;
    FF _FormFac;
    fptype _radius;
    std::vector<Variable *> _AdditionalVars;

  public:
    Lineshape(std::string name,
              Variable *mass,
              Variable *width,
              unsigned int L,
              unsigned int Mpair,
              LS kind                                = LS::BW,
              FF FormFac                             = FF::BL_Prime,
              fptype radius                          = 1.5,
              std::vector<Variable *> AdditionalVars = std::vector<Variable *>());

    bool operator==(const Lineshape &L) const {
        if(_AdditionalVars.size() != L._AdditionalVars.size())
            return false;

        bool addvar = true;

        for(int i = 0; i < _AdditionalVars.size(); ++i) {
            addvar = addvar and (L._AdditionalVars[i]->getValue() == _AdditionalVars[i]->getValue());
        }

        return addvar and (L.getName() == getName() and L._mass->getValue() == _mass->getValue()
                           and L._width->getValue() == _width->getValue()
                           and L._L == _L
                           and L._Mpair == _Mpair
                           and L._kind == _kind
                           and L._FormFac == _FormFac);
    }
    Lineshape(std::string name);

    void setConstantIndex(unsigned int idx) { host_indices[parameters + 1] = idx; }
};

class Amplitude {
    friend class DPPdf;
    friend class TDDP4;

  public:
    Amplitude(std::string uniqueDecayStr,
              Variable *ar,
              Variable *ai,
              std::vector<Lineshape *> LS,
              std::vector<SpinFactor *> SF,
              unsigned int nPerm = 1);
    bool operator==(const Amplitude &A) const;

  private:
    std::string _uniqueDecayStr;
    Variable *_ar;
    Variable *_ai;
    std::vector<SpinFactor *> _SF;
    std::vector<Lineshape *> _LS;
    unsigned int _nPerm;
};
} // namespace GooFit
