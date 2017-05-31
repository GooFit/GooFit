#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class MappedPdf : public GooPdf {
  public:
    MappedPdf(std::string n, GooPdf *m, std::vector<GooPdf *> &t);
    // Map function m must be custom written to correspond to order of function list t.
    __host__ fptype normalize() const override;

  private:
};
} // namespace GooFit
