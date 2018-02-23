#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class NovosibirskPdf : public GooPdf {
  public:
    NovosibirskPdf(std::string n, Observable _x, Variable m, Variable s, Variable t);

    __host__ void recursiveSetIndices() override;

  private:
  private:
};
} // namespace GooFit
