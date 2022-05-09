#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

namespace GooFit {
namespace Lineshapes {

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

} // namespace Lineshapes
} // namespace GooFit
