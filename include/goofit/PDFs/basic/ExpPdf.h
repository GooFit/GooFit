#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ExpPdf : public GooPdf {
  public:
    ExpPdf(std::string n, Observable _x, Variable alpha);
    ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights);
    ExpPdf(std::string n, Observable _x, Variable alpha, Variable offset);
    ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights, Variable offset);
    __host__ fptype integrate(fptype lo, fptype hi) const override;
    __host__ bool hasAnalyticIntegral() const override { return (1 == host_parameters[parametersIdx]); }

    __host__ virtual void recursiveSetIndices();

  private:
    int ExpType;
};
} // namespace GooFit
