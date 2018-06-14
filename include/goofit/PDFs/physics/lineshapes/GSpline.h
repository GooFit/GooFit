#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

#include <tuple>

namespace GooFit {
namespace Lineshapes {

using spline_t = std::tuple<fptype, fptype, unsigned int>;

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

} // namespace GooFit
