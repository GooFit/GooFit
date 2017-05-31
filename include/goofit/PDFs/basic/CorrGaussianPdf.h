#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class CorrGaussianPdf : public GooPdf {
  public:
    CorrGaussianPdf(std::string n,
                    Variable *_x,
                    Variable *_y,
                    Variable *mean1,
                    Variable *sigma1,
                    Variable *mean2,
                    Variable *sigma2,
                    Variable *correlation);

  private:
};

} // namespace GooFit
