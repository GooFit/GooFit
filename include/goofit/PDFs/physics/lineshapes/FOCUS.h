#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

namespace GooFit {
namespace Lineshapes {

class FOCUS final : public Lineshape {
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

    virtual bool isEqualByValue(const Lineshape &other) const override;
};

} // namespace Lineshapes
} // namespace GooFit
