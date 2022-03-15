// This holds anything that does not have explicit CUDA in it.

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/Version.h>

#ifdef ROOT_FOUND
#include <TH1D.h>
#endif

namespace GooFit {

__host__ auto GooPdf::integrate(fptype, fptype) const -> fptype { return 0; }

__host__ auto GooPdf::hasAnalyticIntegral() const -> bool { return false; }

__host__ void GooPdf::setMetrics() {
    logger = std::make_shared<MetricTaker>(this, getMetricPointer(fitControl->getMetric()));
}

/// This collects the number of variables for the thrust call
/// -(n+2) for binned evalute's
__host__ auto GooPdf::get_event_size() const -> int {
    int numVars = observablesList.size();

    if(fitControl->binnedFit())
        return -(numVars + 2);
    else
        return numVars;
}

auto GooPdf::get_bin_grid_size() const -> int {
    size_t totalBins = 1;

    for(const Observable &v : observablesList) {
        totalBins *= integrationBins > 0 ? integrationBins : v.getNumBins();
    }

    return totalBins;
}

auto GooPdf::get_bin_grid_volume() const -> fptype {
    fptype ret = 1.0;

    for(const Observable &v : observablesList) {
        ret *= v.getUpperLimit() - v.getLowerLimit();
    }

    return ret;
}

__host__ void GooPdf::initialize() {
    GOOFIT_DEBUG("{} uses pc.incrementIndex({}, {}, {}, {}, {})",
                 reflex_name_,
                 "1?",
                 parametersList.size(),
                 constantsList.size(),
                 observablesList.size(),
                 "1?");

    if(!fitControl)
        setFitControl(std::make_shared<UnbinnedNllFit>());

    // MetricTaker must be created after PdfBase initialisation is done.
    PdfBase::initializeIndices();
    setMetrics();
}

__host__ void GooPdf::scan(Observable var, std::vector<fptype> &values) {
    fptype step = var.getUpperLimit();
    step -= var.getLowerLimit();
    step /= var.getNumBins();
    values.clear();

    for(fptype v = var.getLowerLimit() + 0.5 * step; v < var.getUpperLimit(); v += step) {
        var.setValue(v);
        copyParams();
        fptype curr = calculateNLL();
        values.push_back(curr);
    }
}

// TODO: is this needed?
__host__ void GooPdf::setParameterConstantness(bool constant) {
    std::vector<Variable> pars = getParameters();

    for(Variable &p : pars) {
        p.setFixed(constant);
    }
}

__host__ auto GooPdf::makeGrid() -> UnbinnedDataSet {
    std::vector<Observable> ret = getObservables();

    UnbinnedDataSet grid{ret};
    grid.fillWithGrid();

    return grid;
}

__host__ void GooPdf::setFitControl_A(std::shared_ptr<FitControl> fc, std::string caller) {
    for(auto &component : components) {
        auto componentName = component->getPdfName();
        component->setFitControl_A(fc, componentName);
    }

    fitControl = fc;

    setMetrics();

    setIndices();
}
__host__ void GooPdf::setFitControl(std::shared_ptr<FitControl> fc) {
    for(auto &component : components) {
        auto componentName = component->getPdfName();
        component->setFitControl(fc);
    }

    fitControl = fc;

    setMetrics();

    setIndices();
}

#ifdef ROOT_FOUND
__host__ auto GooPdf::plotToROOT(Observable var, double normFactor, std::string name) -> TH1D * {
    if(name.empty())
        name = getName() + "_hist";

    auto ret = new TH1D(name.c_str(), "", var.getNumBins(), var.getLowerLimit(), var.getUpperLimit());
    std::vector<fptype> binValues = evaluateAtPoints(var);

    double pdf_int = 0;

    for(int i = 0; i < var.getNumBins(); ++i) {
        pdf_int += binValues[i];
    }

    for(int i = 0; i < var.getNumBins(); ++i)
        ret->SetBinContent(i + 1, binValues[i] * normFactor / pdf_int / var.getBinSize());
    return ret;
}
#endif

} // namespace GooFit
