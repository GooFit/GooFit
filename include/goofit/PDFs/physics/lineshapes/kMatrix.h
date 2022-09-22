#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

#include <array>

#define NPOLES 5
#define NCHANNELS 5

namespace GooFit {
namespace Lineshapes {

class kMatrix : public Lineshape {
  public:
    kMatrix(std::string name,
            Variable a_r,
            Variable a_i,
            Variable sA0,
            Variable sA,
            Variable s0_prod,
            Variable s0_scatt,
            std::vector<Variable> beta_r,
            std::vector<Variable> beta_i,
            std::vector<Variable> f_prod_r,
            std::vector<Variable> f_prod_i,
            std::vector<Variable> fscat,
            std::vector<Variable> poles,
            unsigned int L,
            unsigned int Mpair,
            FF FormFac,
            fptype radius);

    ~kMatrix() override = default;
};

} // namespace Lineshapes
} // namespace GooFit
