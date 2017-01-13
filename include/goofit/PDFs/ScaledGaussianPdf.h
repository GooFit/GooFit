#ifndef SCALEDGAUSSIAN_PDF_HH
#define SCALEDGAUSSIAN_PDF_HH

#include "goofit/PDFs/GooPdf.h"

class ScaledGaussianPdf : public GooPdf {
public:
    ScaledGaussianPdf(std::string n, Variable* _x, Variable* m, Variable* s, Variable* d, Variable* e);
    __host__ virtual bool hasAnalyticIntegral() const {
        return false;
    }



private:

};

#endif
