#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class BWPdf : public GooPdf {
  public:
    BWPdf(std::string n, Observable _x, Variable m, Variable s);
    __host__ void recursiveSetIndices() override;

  private:
};
} // namespace GooFit
