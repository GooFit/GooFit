#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

#include <vector>

namespace GooFit {
namespace Lineshapes {

class GLASS final : public Lineshape {
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

    virtual bool isEqualByValue(const Lineshape &other) const override;

    virtual std::ostream &print(std::ostream &out) const override;
};

} // namespace Lineshapes
} // namespace GooFit
