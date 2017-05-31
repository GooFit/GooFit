#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

// Transforms ND coordinates into a single bin number.
class BinTransformPdf : public GooPdf {
  public:
    BinTransformPdf(std::string n,
                    std::vector<Variable *> obses,
                    std::vector<fptype> limits,
                    std::vector<fptype> binSizes,
                    std::vector<int> numBins);

  private:
};
} // namespace GooFit
