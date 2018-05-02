#include <goofit/Color.h>
#include <goofit/Error.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/Log.h>
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>

#include <algorithm>
#include <set>
#include <utility>

#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
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
fptype host_parameters[maxParams];
fptype host_constants[maxParams];
fptype host_observables[maxParams];
fptype host_normalisations[maxParams];

int host_callnumber     = 0;
int totalParameters     = 0;
int totalConstants      = 0;
int totalObservables    = 0;
int totalNormalisations = 0;

__host__ void PdfBase::checkInitStatus(std::vector<std::string> &unInited) const {
    if(!properlyInitialised)
        unInited.push_back(getName());

    for(auto component : components) {
        component->checkInitStatus(unInited);
    }
}

__host__ void PdfBase::recursiveSetNormalisation(fptype norm) const {
    host_normalisations[normalIdx + 1] = norm;

    for(auto component : components) {
        component->recursiveSetNormalisation(norm);
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
    if(totalConstants + amount >= maxParams)
        throw GooFit::GeneralError(
            "totalConstants {} + amount {} can not be more than {}", totalConstants, amount, maxParams);
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
    setData(data);
    FitManager fitter{this};
    fitter.setVerbosity(verbosity);
    return fitter.fit();
}
} // namespace GooFit
