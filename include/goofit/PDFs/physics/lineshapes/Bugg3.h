#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

namespace GooFit {
namespace Lineshapes {

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

} // namespace Lineshapes
} // namespace GooFit
