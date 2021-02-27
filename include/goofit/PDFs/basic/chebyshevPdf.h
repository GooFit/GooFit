
#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class chebyshevPdf : public GooPdf {
  public:
    chebyshevPdf(std::string n, Observable _x, std::vector<Variable> weights, unsigned int max);
};
__device__ fptype chebyshev(int max, fptype x);
}
