#pragma once

#include <goofit/PDFs/physics/AmpComponent.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/detail/Complex.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

#define NCHANNELS 5

namespace GooFit {
/**
Represents a resonance-shape parametrization, the
\f$B_i\f$ that appear in the equations for GooFit::Amp3Body,
GooFit::Amp3Body_IS, and GooFit::Amp3Body_TD. Canonically a relativistic
Breit-Wigner. The constructor takes the real and imaginary parts of
the coefficient \f$\alpha\f$ (note that this is actually used by the
containing function), and additional parameters depending on which
function the resonance is modeled by:

-   Relativistic Breit-Wigner: Mass, width, spin, and cyclic index.
    The two last are integer constants. Only spins 0, 1, and 2 are
    supported.

-   Gounaris-Sakurai parametrization: Spin, mass, width, and cyclic
    index. Notice that this is the same list as for the relativistic
    BW, just a different order.

-   Nonresonant component (ie, constant across the Dalitz plot):
    Nothing additional.

-   Gaussian: Mean and width of the Gaussian, cyclic index. Notice
    that the Gaussian takes the mass \f$m_{12,13,23}\f$ as its argument,
    not the squared mass \f$m^2_{12,13,23}\f$ like the other
    parametrizations.

Service class intended to hold parametrisations of
resonances on Dalitz plots. Don't try to use this
as a standalone PDF! It should only be used as a
component in one of the friend classes. It extends
GooPdf so as to take advantage of the
infrastructure, but will crash if used on its own.
**/
class ResonancePdf : public AmpComponent {
    friend class Amp3Body;
    friend class Amp3BodySqDP;
    friend class Amp3Body_TD;
    friend class Amp3Body_IS;

  public:
    ~ResonancePdf() override = default;

    __host__ virtual void recalculateCache() const {}

    __host__ auto get_amp_real() const -> Variable { return amp_real; }
    __host__ auto get_amp_img() const -> Variable { return amp_imag; }

  protected:
    /// Special constructor that subclasses use
    ResonancePdf(std::string pdf_name, std::string name, Variable ar, Variable ai)
        : AmpComponent("Resonances::" + pdf_name, name)
        , amp_real(ar)
        , amp_imag(ai) {}

    Variable amp_real;
    Variable amp_imag;

    std::vector<unsigned int> pindices;

    std::vector<fptype> host_constants;
};

} // namespace GooFit
