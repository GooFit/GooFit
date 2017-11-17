#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class NovosibirskPdf : public GooPdf {
  public:
    NovosibirskPdf(std::string n, Observable _x, Variable m, Variable s, Variable t);

  private:
};
} // namespace GooFit
