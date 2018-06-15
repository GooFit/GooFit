#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/Amp3Body_TD.h>

namespace GooFit {

/**
Contains a cyclic
index (either `PAIR_12`, `PAIR_13`, or `PAIR_23`) and the lower and
upper bounds of the veto region.
**/

struct VetoInfo {
    Variable minimum;
    Variable maximum;
    DaughterPair cyclic_index;
    VetoInfo(Variable minimum, Variable maximum, DaughterPair cyclic_index)
        : minimum(minimum)
        , maximum(maximum)
        , cyclic_index(cyclic_index) {}
};

/**
Tests whether a point is in a particular region of
the Dalitz plot, and returns zero if so, one otherwise. Intended for
use as part of an efficiency function, excluding particular
regions - canonically the one containing the \f$K^0\to\pi\pi\f$ decay,
as a large source of backgrounds that proved hard to model. The
constructor takes the squared-mass variables \f$m_{12}\f$ and \f$m_{13}\f$,
the masses (contained in GooFit::Variables) of the mother and three
daughter particles involved in the decay, and a `vector` of
GooFit::VetoInfo objects.

The GooFit::VetoInfo objects just contain a cyclic
index (either `PAIR_12`, `PAIR_13`, or `PAIR_23`) and the lower and
upper bounds of the veto region.
**/

class DalitzVetoPdf : public GooPdf {
  public:
    __host__ DalitzVetoPdf(std::string n,
                           Observable _x,
                           Observable _y,
                           Variable motherM,
                           Variable d1m,
                           Variable d2m,
                           Variable d3m,
                           std::vector<VetoInfo> vetos);
};
} // namespace GooFit
