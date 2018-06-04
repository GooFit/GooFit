#include <goofit/Color.h>
#include <goofit/Error.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/Log.h>
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>

#include <algorithm>
#include <random>
#include <set>
#include <utility>

#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/detail/CompilerFeatures.h>

#include <Minuit2/FunctionMinimum.h>

namespace {

template <typename T>
bool find_in(std::vector<T> list, T item) {
    return std::find_if(std::begin(list), std::end(list), [item](T p) { return p == item; }) != std::end(list);
}
} // namespace
namespace GooFit {

fptype *dev_event_array;
fptype host_parameters[GOOFIT_MAXPAR];
fptype host_constants[GOOFIT_MAXPAR];
fptype host_observables[GOOFIT_MAXPAR];
fptype host_normalizations[GOOFIT_MAXPAR];

int host_callnumber     = 0;
int totalParameters     = 0;
int totalConstants      = 0;
int totalObservables    = 0;
int totalNormalizations = 0;

__host__ void PdfBase::checkInitStatus(std::vector<std::string> &unInited) const {
    if(!properlyInitialised)
        unInited.push_back(getName());

    for(auto component : components) {
        component->checkInitStatus(unInited);
    }
}

__host__ void PdfBase::recursiveSetNormalization(fptype norm) const {
    host_normalizations[normalIdx + 1] = norm;

    for(auto component : components) {
        component->recursiveSetNormalization(norm);
    }
}

__host__ unsigned int PdfBase::registerParameter(Variable var) {
    // if(find_in(parametersList, var))
    //    throw GeneralError("This var {} was registered twice", var.getName());

    parametersList.push_back(var);

    return parametersList.size() - 1; // TODO: Make void
}

__host__ unsigned int PdfBase::registerConstant(fptype value) {
    constantsList.push_back(value);
    // Return the index that this variable is stored at, not the total!
    return constantsList.size() - 1;
}

__host__ void PdfBase::unregisterParameter(Variable var) {
    GOOFIT_DEBUG("{}: Removing {}", getName(), var.getName());

    auto pos = std::find(parametersList.begin(), parametersList.end(), var);

    if(pos != parametersList.end())
        parametersList.erase(pos);

    for(PdfBase *comp : components) {
        comp->unregisterParameter(var);
    }
}

__host__ std::vector<Variable> PdfBase::getParameters() const {
    std::vector<Variable> ret;
    for(const Variable &param : parametersList)
        ret.push_back(param);

    for(const PdfBase *comp : components) {
        for(const Variable &sub_comp : comp->getParameters())
            if(!find_in(ret, sub_comp))
                ret.push_back(sub_comp);
    }

    return ret;
}

__host__ Variable *PdfBase::getParameterByName(std::string n) {
    for(Variable &p : parametersList) {
        if(p.getName() == n)
            return &p;
    }

    for(auto component : components) {
        Variable *cand = component->getParameterByName(n);

        if(cand != nullptr)
            return cand;
    }

    return nullptr;
}

__host__ std::vector<Observable> PdfBase::getObservables() const {
    std::vector<Observable> ret;
    for(const Observable &obs : observablesList)
        ret.push_back(obs);

    for(const PdfBase *comp : components) {
        for(const Observable &sub_comp : comp->getObservables())
            if(!find_in(ret, sub_comp))
                ret.push_back(sub_comp);
    }

    return ret;
}

GOOFIT_DEPRECATED __host__ unsigned int PdfBase::registerConstants(unsigned int amount) {
    if(totalConstants + amount >= GOOFIT_MAXPAR)
        throw GooFit::GeneralError(
            "totalConstants {} + amount {} can not be more than {}", totalConstants, amount, GOOFIT_MAXPAR);
    cIndex = totalConstants;
    totalConstants += amount;
    return cIndex;
}

void PdfBase::registerObservable(Observable obs) {
    if(find_in(observablesList, obs))
        return;

    observablesList.push_back(obs);
}

void PdfBase::setupObservables() {
    int counter = 0;
    for(auto &i : observablesList) {
        i.setIndex(counter);
        counter++;
    }

    GOOFIT_DEBUG("SetupObservables {} ", counter);
}

__host__ void PdfBase::setIntegrationFineness(int i) {
    integrationBins = i;
    generateNormRange();
}

__host__ bool PdfBase::parametersChanged() const {
    return std::any_of(
        std::begin(parametersList), std::end(parametersList), [](const Variable &v) { return v.getChanged(); });
}

__host__ void PdfBase::setNumPerTask(PdfBase *p, const int &c) {
    if(p == nullptr)
        return;

    m_iEventsPerTask = c;
}

__host__ ROOT::Minuit2::FunctionMinimum PdfBase::fitTo(DataSet *data, int verbosity) {
    auto old = getData();
    setData(data);
    auto funmin = fit(verbosity);
    setData(old);
    return funmin;
}

__host__ ROOT::Minuit2::FunctionMinimum PdfBase::fit(int verbosity) {
    FitManager fitter{this};
    fitter.setVerbosity(verbosity);
    return fitter.fit();
}

void PdfBase::fillMCDataSimple(size_t events, unsigned int seed) {
    // Setup bins
    if(observablesList.size() != 1)
        throw GeneralError("You can only fill MC on a 1D dataset with a simple fill");

    Observable var = observablesList.at(0);

    UnbinnedDataSet data{var};
    data.fillWithGrid();
    auto origdata = getData();

    if(origdata == nullptr)
        throw GeneralError("Can't run on a PDF with no DataSet to fill!");

    setData(&data);
    std::vector<double> pdfValues = dynamic_cast<GooPdf *>(this)->getCompProbsAtDataPoints()[0];

    // Setup random numbers
    if(seed == 0) {
        std::random_device rd;
        seed = rd();
    }
    std::mt19937 gen(seed);

    // Uniform distribution
    std::uniform_real_distribution<> unihalf(-.5, .5);
    std::uniform_real_distribution<> uniwhole(0.0, 1.0);

    // CumSum in other languages
    std::vector<double> integral(pdfValues.size());
    std::partial_sum(pdfValues.begin(), pdfValues.end(), integral.begin());

    // Make this a 0-1 fraction by dividing by the end value
    std::for_each(integral.begin(), integral.end(), [&integral](double &val) { val /= integral.back(); });

    for(size_t i = 0; i < events; i++) {
        double r = uniwhole(gen);

        // Binary search for integral[cell-1] < r < integral[cell]
        size_t j = std::lower_bound(integral.begin(), integral.end(), r) - integral.begin();

        // Fill in the grid randomly
        double varValue = data.getValue(var, j) + var.getBinSize() * unihalf(gen);

        var.setValue(varValue);
        origdata->addEvent();
    }

    setData(origdata);
}

std::ostream &operator<<(std::ostream &out, const PdfBase &pdf) {
    out << "GooPdf(\"" << pdf.getName() << "\") :\n";
    out << "  Device function: " << pdf.reflex_name_ << "\n";

    out << "  Observable" << (pdf.observablesList.size() > 1 ? "s" : "") << ": ";
    for(const auto &item : pdf.observablesList)
        out << item.getName() << (item == pdf.observablesList.back() ? "" : ", ");
    out << "\n";

    out << "  Parameter" << (pdf.parametersList.size() > 1 ? "s" : "") << ": ";
    for(const auto &item : pdf.parametersList)
        out << item.getName() << (item == pdf.parametersList.back() ? "" : ", ");
    out << "\n";

    return out;
}

} // namespace GooFit
