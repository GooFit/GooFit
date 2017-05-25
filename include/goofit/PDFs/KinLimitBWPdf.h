#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {


class KinLimitBWPdf : public GooPdf {

public:
    KinLimitBWPdf(std::string n, Variable* _x, Variable* m, Variable* s);
    __host__ virtual bool hasAnalyticIntegral() const {
        return false;
    }
    __host__ void setMasses(fptype bigM, fptype smallM);

private:

};
} // namespace GooFit

