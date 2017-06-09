#pragma once

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/physics/TddpPdf.h"

namespace GooFit {

struct VetoInfo {
    DaughterPair cyclic_index;
    Variable *minimum;
    Variable *maximum;
};

class DalitzVetoPdf : public GooPdf {
  public:
    __host__ DalitzVetoPdf(std::string n,
                           Variable *_x,
                           Variable *_y,
                           Variable *motherM,
                           Variable *d1m,
                           Variable *d2m,
                           Variable *d3m,
                           std::vector<VetoInfo *> vetos);

  private:
};
} // namespace GooFit
