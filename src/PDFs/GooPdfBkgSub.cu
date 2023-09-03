#include <goofit/PDFs/GooPdfBkgSub.h>

namespace GooFit {

GooPdfBkgSub::GooPdfBkgSub(GooPdf basePdf, UnbinnedDataSet *bkgEvents, Variable nBkg)
    : GooPdf(basePdf)
    , m_bkgEvents(bkgEvents) {
    registerParameter(nBkg);
}

__host__ double GooPdfBkgSub::calculateNLL() {
    // Calculate the NLL for the fit data.
    double ret = GooPdf::calculateNLL();

    // Subtract the background contribution to the NLL.
    auto oldData = GooPdf::getData();
    GooPdf::setData(m_bkgEvents);
    size_t nBkgEvents = m_bkgEvents->numberOfEvents();
    fptype nBkgFit    = parametersList.back().getValue();
    ret -= GooPdf::calculateNLL() * nBkgFit / nBkgEvents;

    // Reset the fit data.
    GooPdf::setData(oldData);

    return ret;
}

} // namespace GooFit
