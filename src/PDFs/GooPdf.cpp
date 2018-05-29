// This holds anything that does not have explicit CUDA in it.

#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/Version.h>

#ifdef ROOT_FOUND
#include <TH1D.h>
#endif

namespace GooFit {

void *host_function_table[GOOFIT_MAXFUNC];
void *host_fcn_ptr                = nullptr;
unsigned int num_device_functions = 0;

__host__ fptype GooPdf::integrate(fptype, fptype) const { return 0; }

__host__ bool GooPdf::hasAnalyticIntegral() const { return false; }

__host__ void GooPdf::setMetrics() {
    logger = std::make_shared<MetricTaker>(this, getMetricPointer(fitControl->getMetric()));
}

__host__ void GooPdf::initialize() {
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

__host__ UnbinnedDataSet GooPdf::makeGrid() {
    std::vector<Observable> ret = getObservables();

    UnbinnedDataSet grid{ret};
    grid.fillWithGrid();

    return grid;
}

__host__ void GooPdf::setFitControl(std::shared_ptr<FitControl> fc) {
    for(auto &component : components) {
        component->setFitControl(fc);
    }

    fitControl = fc;

    setMetrics();

    setIndices();
}

#ifdef ROOT_FOUND
__host__ TH1D *GooPdf::plotToROOT(Observable var, double normFactor, std::string name) {
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
