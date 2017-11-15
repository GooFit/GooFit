#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class VoigtianPdf : public GooPdf {
  public:
    VoigtianPdf(std::string n, Observable _x, Variable m, Variable s, Variable w);

  private:
};
} // namespace GooFit
