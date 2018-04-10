/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/Variable.h>
#include <goofit/Version.h>

#include <array>

#include <goofit/detail/Complex.h>

namespace GooFit {

class SpinFactor;

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

    std::vector<unsigned int> pindices{0};

    /// Protected constructor, only for subclasses to use
    Lineshape(std::string name);

  public:
    ~Lineshape() override = default;

    void setConstantIndex(unsigned int idx) { host_indices[parameters + 1] = idx; }

    bool operator==(const Lineshape &L) const { return (L.getName() == getName()); }
};

namespace Lineshapes {

using spline_t = std::tuple<fptype, fptype, unsigned int>;

class RBW : public Lineshape {
  public:
    RBW(std::string name,
        Variable mass,
        Variable width,
        unsigned int L,
        unsigned int Mpair,
        FF FormFac    = FF::BL_Prime,
        fptype radius = 1.5);

    ~RBW() override = default;
};

class One : public Lineshape {
  public:
    One(std::string name,
        Variable mass,
        Variable width,
        unsigned int L,
        unsigned int Mpair,
        FF FormFac    = FF::BL_Prime,
        fptype radius = 1.5);

    ~One() override = default;
};

class LASS : public Lineshape {
  public:
    LASS(std::string name,
         Variable mass,
         Variable width,
         unsigned int L,
         unsigned int Mpair,
         FF FormFac    = FF::BL_Prime,
         fptype radius = 1.5);

    ~LASS() override = default;
};

class NonRes : public Lineshape {
  public:
    NonRes(std::string name,
           Variable mass,
           Variable width,
           unsigned int L,
           unsigned int Mpair,
           FF FormFac    = FF::BL_Prime,
           fptype radius = 1.5);

    ~NonRes() override = default;
};

class Bugg : public Lineshape {
  public:
    Bugg(std::string name,
         Variable mass,
         Variable width,
         unsigned int L,
         unsigned int Mpair,
         FF FormFac    = FF::BL_Prime,
         fptype radius = 1.5);

    ~Bugg() override = default;
};

class Bugg3 : public Lineshape {
  public:
    Bugg3(std::string name,
          Variable mass,
          Variable width,
          unsigned int L,
          unsigned int Mpair,
          FF FormFac    = FF::BL_Prime,
          fptype radius = 1.5);

    ~Bugg3() override = default;
};

class Flatte : public Lineshape {
  public:
    Flatte(std::string name,
           Variable mass,
           Variable width,
           unsigned int L,
           unsigned int Mpair,
           FF FormFac    = FF::BL_Prime,
           fptype radius = 1.5);

    ~Flatte() override = default;
};

class SBW : public Lineshape {
  public:
    SBW(std::string name,
        Variable mass,
        Variable width,
        unsigned int L,
        unsigned int Mpair,
        FF FormFac    = FF::BL_Prime,
        fptype radius = 1.5);

    ~SBW() override = default;
};

class FOCUS : public Lineshape {
  public:
    enum class Mod { Kpi = 0, KEta, I32 };
    FOCUS(std::string name,
          Mod mod,
          Variable mass,
          Variable width,
          unsigned int L,
          unsigned int Mpair,
          FF FormFac    = FF::BL_Prime,
          fptype radius = 1.5);

    ~FOCUS() override = default;
};

#if GOOFIT_KMATRIX
class kMatrix : public Lineshape {
  public:
    kMatrix(std::string name,
            unsigned int pterm, //< 0 or 1
            bool is_pole,       //< False for prod
            Variable sA0,
            Variable sA,
            Variable s0_prod,
            Variable s0_scatt,
            std::array<Variable, 5> f,
            std::array<Variable, 5 * 6> poles,
            Variable mass,
            Variable width,
            unsigned int L,
            unsigned int Mpair,
            FF FormFac    = FF::BL_Prime,
            fptype radius = 1.5);

    ~kMatrix() override = default;
};
#endif

class GLASS : public Lineshape {
  public:
    GLASS(std::string name,
          Variable mass,
          Variable width,
          unsigned int L,
          unsigned int Mpair,
          FF FormFac,
          fptype radius,
          std::vector<Variable> AdditionalVars);

    ~GLASS() override = default;
};

/// A spline implementaiton for the width (Gamma = G)
class GSpline : public Lineshape {
  public:
    /// The spline
    GSpline(std::string name,
            Variable mass,
            Variable width,
            unsigned int L,
            unsigned int Mpair,
            FF FormFac,
            fptype radius,
            std::vector<Variable> AdditionalVars,
            spline_t SplineInfo);

    ~GSpline() override = default;
};
} // namespace Lineshapes

class Amplitude {
    friend class DPPdf;
    friend class TDDP4;

  public:
    Amplitude(std::string uniqueDecayStr,
              Variable ar,
              Variable ai,
              std::vector<Lineshape *> LS,
              std::vector<SpinFactor *> SF,
              unsigned int nPerm = 1);

    bool operator==(const Amplitude &A) const;

  private:
    std::string _uniqueDecayStr;
    Variable _ar;
    Variable _ai;
    std::vector<SpinFactor *> _SF;
    std::vector<Lineshape *> _LS;
    unsigned int _nPerm;
};
} // namespace GooFit
