#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <vector>

namespace GooFit {

// Transforms ND coordinates into a single bin number.
class VariableBinTransform1DPdf : public GooPdf {
  public:
    VariableBinTransform1DPdf(std::string n, Observable _x, std::vector<fptype> binlimits);
};

} // namespace GooFit
