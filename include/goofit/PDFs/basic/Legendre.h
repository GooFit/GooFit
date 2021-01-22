#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class LegendrePdf : public GooPdf {
  public:
    LegendrePdf(std::string n, Observable _x, std::vector<Variable> weights, unsigned int max);
};
__device__ fptype Legendre(int max, fptype x);
__device__ int Factorial(int n);
__device__ fptype BinomCoeff(int n, int k);
} // namespace GooFit
