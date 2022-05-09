
#pragma once

#include <goofit/PDFs/physics/AmpComponent.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/Variable.h>
#include <goofit/Version.h>
#include <goofit/detail/Complex.h>

#include <array>

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

class Lineshape : public AmpComponent {
  protected:
    friend class Amp4Body;
    friend class Amp4Body_TD;

    /// Protected constructor, only for subclasses to use
    Lineshape(std::string pdf_name, std::string name, unsigned int L, unsigned int Mpair, FF FormFac, fptype radius);

    unsigned int _L;
    unsigned int _Mpair;
    FF _FormFac;
    fptype _radius;

  public:
    ~Lineshape() override = default;

    auto operator==(const Lineshape &L) const -> bool { return L.getName() == getName(); }
};

} // namespace GooFit
