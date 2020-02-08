#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <array>

#define NPOLES 5
#define NCHANNELS 5

namespace GooFit {
namespace Resonances {

class kMatrix : public ResonancePdf {
  public:
    kMatrix(std::string name,
            unsigned int pterm, //< 0 or 1
            bool is_pole,       //< False for prod
            Variable a_r,
            Variable a_i,
            Variable sA0,
            Variable sA,
            Variable s0_prod,
            Variable s0_scatt,
            std::vector<Variable> fscat,
            std::vector<Variable> poles,
            unsigned int L,
            unsigned int Mpair);

    ~kMatrix() override = default;
};

} // namespace Resonances
} // namespace GooFit
