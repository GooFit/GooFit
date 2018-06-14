#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

#include <vector>

namespace GooFit {
namespace Lineshapes {

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

} // namespace Lineshapes
} // namespace GooFit
