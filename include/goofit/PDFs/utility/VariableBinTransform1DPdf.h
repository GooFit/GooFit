#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <vector>

namespace GooFit {

// Transforms ND coordinates into a single bin number.
class VariableBinTransform1DPdf : public GooPdf {
  public:
    VariableBinTransform1DPdf(std::string n, Variable *_x, std::vector<fptype> binlimits);
};

} // GooFit
