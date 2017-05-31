#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class TrigThresholdPdf : public GooPdf {
  public:
    TrigThresholdPdf(
        std::string n, Variable *_x, Variable *thresh, Variable *trigConst, Variable *linConst, bool upper = true);
    TrigThresholdPdf(std::string n,
                     Variable *_x,
                     Variable *_y,
                     Variable *thresh,
                     Variable *trigConst,
                     Variable *linConst,
                     Variable *massConstant,
                     bool upper);

  private:
};
} // namespace GooFit
