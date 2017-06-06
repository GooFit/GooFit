#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {


class GaussianPdf : public GooPdf {
public:
    GaussianPdf(std::string n, Variable* _x, Variable* m, Variable* s);
    __host__ fptype integrate(fptype lo, fptype hi) const;
    __host__ virtual bool hasAnalyticIntegral() const {
        return true;
    }

protected:
    __host__ virtual void recursiveSetIndices ();

private:

};
} // namespace GooFit

