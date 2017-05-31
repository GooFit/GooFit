#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {

class PolynomialPdf : public GooPdf {
  public:
    PolynomialPdf(std::string n,
                  Variable *_x,
                  std::vector<Variable *> weights,
                  Variable *x0              = nullptr,
                  unsigned int lowestDegree = 0);
    PolynomialPdf(std::string n,
                  std::vector<Variable *> obses,
                  std::vector<Variable *> coeffs,
                  std::vector<Variable *> offsets,
                  unsigned int maxDegree);
    __host__ fptype integrate(fptype lo, fptype hi) const override;
    //__host__ virtual bool hasAnalyticIntegral () const {return (1 == observables.size());}
    __host__ fptype getCoefficient(int coef) const;

  private:
    Variable *center;
};
} // namespace GooFit
