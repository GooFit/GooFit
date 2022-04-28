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

#include <Minuit2/FunctionMinimum.h>

namespace {

template <typename T>
auto find_in(std::vector<T> list, T item) -> bool {
    return std::find_if(std::begin(list), std::end(list), [item](T p) { return p == item; }) != std::end(list);
}
} // namespace

namespace GooFit {

__host__ void PdfBase::recursiveSetNormalization(fptype norm, bool subpdf) {
    host_normalizations.at(normalIdx + 1) = norm;
    cachedNormalization                   = norm;

    for(auto component : components) {
        component->recursiveSetNormalization(norm, true);
    }
    if(!subpdf)
        host_normalizations.sync(d_normalizations);
}

__host__ void PdfBase::registerParameter(Variable var) { parametersList.push_back(var); }

__host__ void PdfBase::registerConstant(fptype value) { constantsList.push_back(value); }

__host__ void PdfBase::unregisterParameter(Variable var) {
    GOOFIT_DEBUG("{}: Removing {}", getName(), var.getName());

    auto pos = std::find(parametersList.begin(), parametersList.end(), var);

    if(pos != parametersList.end())
        parametersList.erase(pos);

    for(PdfBase *comp : components) {
        comp->unregisterParameter(var);
    }
}

__host__ auto PdfBase::getParameters() const -> std::vector<Variable> {
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

__host__ auto PdfBase::getParameterByName(std::string n) -> Variable * {
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

__host__ auto PdfBase::getObservables() const -> std::vector<Observable> {
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

__host__ auto PdfBase::parametersChanged() const -> bool {
    return std::any_of(
        std::begin(parametersList), std::end(parametersList), [](const Variable &v) { return v.getChanged(); });
}

__host__ void PdfBase::setNumPerTask(PdfBase *p, const int &c) {
    if(p == nullptr)
        return;

    m_iEventsPerTask = c;

    // we need to set all children components
    for(auto &component : components)
        component->setNumPerTask(component, c);
}

__host__ auto PdfBase::fitTo(DataSet *data, int verbosity) -> ROOT::Minuit2::FunctionMinimum {
    auto old = getData();
    setData(data);
    auto funmin = fit(verbosity);
    setData(old);
    return funmin;
}

__host__ auto PdfBase::fit(int verbosity) -> ROOT::Minuit2::FunctionMinimum {
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

auto operator<<(std::ostream &out, const PdfBase &pdf) -> std::ostream & {
    out << "GooPdf::" << pdf.getPdfName() << "(\"" << pdf.getName() << "\") :\n";
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

bool PdfBase::areParamsandConstantsEqualByVal(const PdfBase &other) const {
    if(this->getParameters().size() != other.getParameters().size()) {
        return false;
    }
    for(int p = 0; p < this->getParameters().size(); p++) {
        Variable thisParam  = this->getParameters()[p];
        Variable otherParam = other.getParameters()[p];
        if(!thisParam.isEqualNameValLimits(otherParam)) {
            return false;
        }
    }

    if(this->constantsList.size() != other.constantsList.size()) {
        return false;
    }
    for(int c = 0; c < this->constantsList.size(); c++) {
        fptype thisConst  = this->constantsList[c];
        fptype otherConst = other.constantsList[c];
        if(thisConst != otherConst) {
            return false;
        }
    }

    return true;
}

} // namespace GooFit
