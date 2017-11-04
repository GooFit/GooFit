#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class LandauPdf : public GooPdf {
  public:
    LandauPdf(std::string n, Variable *_x, Variable *mpv, Variable *sigma);

    __host__ virtual void recursiveSetIndices();

  private:
  private:
};
} // namespace GooFit
