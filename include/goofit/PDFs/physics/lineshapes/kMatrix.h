#pragma once

#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

#include <array>

#define NPOLES 5
#define NCHANNELS 5

namespace GooFit {
namespace Lineshapes {

class kMatrixNew : public Lineshape {
  public:
    kMatrixNew(std::string name,
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

    ~kMatrixNew() override = default;
};

} // namespace Lineshapes
} // namespace GooFit
