#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ExpGausPdf : public GooPdf {
  public:
    ExpGausPdf(std::string n, Observable _x, Variable m, Variable s, Variable t);

  private:
};
} // namespace GooFit
