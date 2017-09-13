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

#include <array>

#include "goofit/detail/Complex.h"

namespace GooFit {

class SpinFactor;

enum class LS { ONE, BW, Lass, Lass_M3, nonRes, Bugg, Bugg3, Flatte, SBW, Spline, kMatrix};
    
// PDG notation for FF
enum class FF : unsigned int { One = 0, BL, BL_Prime, BL2 };

/// Service class intended to hold parametrisations of
/// resonances on Dalitz plots. Don't try to use this
/// as a standalone PDF! It should only be used as a
/// component in one of the friend classes. It extends
/// GooPdf so as to take advantage of the
/// infrastructure, but will crash if used on its own.
    
class Lineshape : public GooPdf {
protected:
    friend class DPPdf;
    friend class TDDP4;
    
    Variable *_mass;
    Variable *_width;
    unsigned int _L;
    unsigned int _Mpair;
    
    LS _kind;
    FF _FormFac;
    
    fptype _radius;
    
    std::vector<unsigned int> pindices {0};
    
    /// Protected constructor, only for subclasses to use
    Lineshape(
        Variable *,
        std::string name,
        Variable *mass,
        Variable *width,
        unsigned int L,
        unsigned int Mpair,
        LS kind,
        FF FormFac,
        fptype radius);
    
public:
    /// Construct standard versions
    Lineshape(std::string name,
                  Variable *mass,
                  Variable *width,
                  unsigned int L,
                  unsigned int Mpair,
                  LS kind,
                  FF FormFac = FF::BL_Prime,
                  fptype radius = 1.5);

    
    virtual ~Lineshape() = default;
    
    void setConstantIndex(unsigned int idx) { host_indices[parameters + 1] = idx; }
    
    bool operator==(const Lineshape &L) const {
        return  (L.getName() == getName() and L._mass->getValue() == _mass->getValue()
                           and L._width->getValue() == _width->getValue()
                           and L._L == _L
                           and L._Mpair == _Mpair
                           and L._kind == _kind
                           and L._FormFac == _FormFac);
    }

};
    
namespace Lineshapes {
using spline_t = std::tuple<fptype, fptype, unsigned int>;
    
class RBW : public Lineshape {
    public:
    
    RBW(std::string name,
              Variable *mass,
              Variable *width,
              unsigned int L,
              unsigned int Mpair,
              FF FormFac                             = FF::BL_Prime,
              fptype radius                          = 1.5);
    
    virtual ~RBW() = default;
        
};
    
class kMatrix : public Lineshape {
public:
    enum class Channels : unsigned int {
        pipi=0, KK, pipipipi, EtaEta, EtapEta
    };

protected:
    unsigned int pterm;
    Channels channels;
    
public:
    kMatrix(std::string name,
        unsigned int pterm,
        Channels chan,
        Variable* sA0,
        Variable* sA,
        Variable* s0_prod,
        Variable* s0_scatt,
        std::array<Variable*, 5> f,
        std::array<Variable*, 5*6> poles,
        Variable *mass,
        Variable *width,
        unsigned int L,
        unsigned int Mpair,
        FF FormFac                             = FF::BL_Prime,
        fptype radius                          = 1.5);
    
    virtual ~kMatrix() = default;
    
};
    
class LASS : public Lineshape {
    
    protected:
        
        
        std::vector<Variable *> _AdditionalVars;
        
    public:
        
        LASS(std::string name,
                Variable *mass,
                Variable *width,
                unsigned int L,
                unsigned int Mpair,
                FF FormFac,
                fptype radius,
                std::vector<Variable *> AdditionalVars);
        
        virtual ~LASS() = default;
        
    };

class GSpline : public Lineshape {

protected:

    std::vector<Variable *> _AdditionalVars;
    std::vector<Variable *> _Curvature;
    
    /// Min, Max, N for parametrization
    spline_t _SplineInfo;

public:
    
    GSpline(std::string name,
              Variable *mass,
              Variable *width,
              unsigned int L,
              unsigned int Mpair,
              FF FormFac,
              fptype radius,
              std::vector<Variable *> AdditionalVars,
              std::vector<Variable *> Curvatures,
              spline_t SplineInfo);

    virtual ~GSpline() = default;
    
};
    
}

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
