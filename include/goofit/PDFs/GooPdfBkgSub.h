#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class GooPdfBkgSub : public GooPdf {
    /**
    This class implements a version of GooPdf that includes a sample of
    background events that are used to statistically subtract the background
    contributions from the resulting likelihood. This can be used to handle
    backgrounds without a convenient parameterization, such as combinatoric
    background in amplitude analyses.
    **/

  public:
    GooPdfBkgSub(GooPdf basePdf, UnbinnedDataSet *bkgEvents, Variable bkgFraction);

  protected:
    __host__ double calculateNLL() override;

  private:
    UnbinnedDataSet *m_bkgEvents;
};

} // namespace GooFit
