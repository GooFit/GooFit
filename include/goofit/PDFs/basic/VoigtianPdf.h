#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class VoigtianPdf : public GooPdf {
  public:
    VoigtianPdf(std::string n, Variable *_x, Variable *m, Variable *s, Variable *w);

    __host__ virtual void recursiveSetIndices ();

  private:
};
} // namespace GooFit
