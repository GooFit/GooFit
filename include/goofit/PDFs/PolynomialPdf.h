#ifndef POLYNOMIAL_PDF_HH
#define POLYNOMIAL_PDF_HH

#include "goofit/PDFs/GooPdf.h"

class PolynomialPdf : public GooPdf {
public:
    PolynomialPdf(std::string n, Variable* _x, std::vector<Variable*> weights, Variable* x0 = 0,
                  unsigned int lowestDegree = 0);
    PolynomialPdf(string n, vector<Variable*> obses, vector<Variable*> coeffs, vector<Variable*> offsets,
                  unsigned int maxDegree);
    __host__ fptype integrate(fptype lo, fptype hi) const;
    //__host__ virtual bool hasAnalyticIntegral () const {return (1 == observables.size());}
    __host__ fptype getCoefficient(int coef) const;


private:
    Variable* center;
};

#endif
