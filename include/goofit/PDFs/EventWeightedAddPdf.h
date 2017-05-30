#pragma once

#include "goofit/PDFs/GooPdf.h"

namespace GooFit {


// This class is just like AddPdf except that the
// event weights are properties of each event, not variables
// in the fit.
class EventWeightedAddPdf : public GooPdf {
public:

    EventWeightedAddPdf(std::string n, std::vector<Variable*> weights, std::vector<PdfBase*> comps);
    __host__ virtual fptype normalize() const;
    __host__ virtual bool hasAnalyticIntegral() const {
        return false;
    }

protected:

private:

};
} // namespace GooFit

