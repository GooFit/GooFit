#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class Component : public GooPdf {
  public:
    using GooPdf::GooPdf;
};

} // namespace GooFit
