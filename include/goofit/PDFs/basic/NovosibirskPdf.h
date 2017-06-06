#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class NovosibirskPdf : public GooPdf {
public:
    NovosibirskPdf(std::string n, Variable* _x, Variable* m, Variable* s, Variable* t);

    __host__ virtual void recursiveSetIndices ();
private:

  private:
};
} // namespace GooFit
