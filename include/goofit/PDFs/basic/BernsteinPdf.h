#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
A Bernstein polynomial, named after Sergei Natanovich Bernstein, is a polynomial
in the Bernstein form, that is a linear combination of Bernstein basis polynomials.

**/
class BernsteinPdf : public GooPdf {
  public:
    BernsteinPdf(std::string n, Observable vars, std::vector<Variable> coeffs, unsigned int mindeg = 0);
    __host__ auto hasAnalyticIntegral() const -> bool override { return false; }

  private:
    fptype *host_constants;
};
} // namespace GooFit
