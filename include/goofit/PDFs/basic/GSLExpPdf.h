#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class GSLExpPdf : public GooPdf {
  public:
    GSLExpPdf(std::string n, Observable _x, Variable alpha);
};

} // namespace GooFit
