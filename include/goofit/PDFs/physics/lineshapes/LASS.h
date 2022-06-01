#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

namespace GooFit {
namespace Lineshapes {

class LASS final : public Lineshape {
  public:
    LASS(std::string name,
         Variable mass,
         Variable width,
         unsigned int L,
         unsigned int Mpair,
         FF FormFac    = FF::BL_Prime,
         fptype radius = 1.5);

    ~LASS() override = default;

    virtual bool isEqualByValue(const Lineshape &other) const override;
};

} // namespace Lineshapes
} // namespace GooFit
