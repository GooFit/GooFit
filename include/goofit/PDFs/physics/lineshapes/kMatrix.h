#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

#include <array>

namespace GooFit {
namespace Lineshapes {

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

} // namespace Lineshapes
} // namespace GooFit
