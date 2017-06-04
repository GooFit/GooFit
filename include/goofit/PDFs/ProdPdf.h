#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {


class ProdPdf : public GooPdf {
public:

    ProdPdf(std::string n, std::vector<PdfBase*> comps);
    __host__ virtual fptype normalize() const;
    __host__ virtual bool hasAnalyticIntegral() const {
        return false;
    }

protected:
    __host__ virtual void recursiveSetIndices ();

private:
    bool varOverlaps; // True if any components share an observable.
};
} // namespace GooFit

