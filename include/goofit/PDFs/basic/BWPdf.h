#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class BWPdf : public GooPdf {
  public:
    BWPdf(std::string n, Observable _x, Variable m, Variable s);

  private:
};
} // namespace GooFit
