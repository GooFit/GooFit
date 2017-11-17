#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class TrigThresholdPdf : public GooPdf {
  public:
    TrigThresholdPdf(
        std::string n, Observable _x, Variable thresh, Variable trigConst, Variable linConst, bool upper = true);
    TrigThresholdPdf(std::string n,
                     Observable _x,
                     Observable _y,
                     Variable thresh,
                     Variable trigConst,
                     Variable linConst,
                     Variable massConstant,
                     bool upper);

  private:
};
} // namespace GooFit
