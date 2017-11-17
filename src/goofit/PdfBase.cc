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

#include <Minuit2/FunctionMinimum.h>

namespace {

template <typename T>
bool find_in(std::vector<T> list, T item) {
    return std::find_if(std::begin(list), std::end(list), [item](T p) { return p == item; }) != std::end(list);
}
} // namespace

namespace GooFit {

fptype *dev_event_array;
fptype host_normalisation[maxParams];
fptype host_params[maxParams];
unsigned int host_indices[maxParams];

int host_callnumber = 0;
int totalParams     = 0;
int totalConstants  = 1; // First constant is reserved for number of events.

__host__ void PdfBase::checkInitStatus(std::vector<std::string> &unInited) const {
    if(!properlyInitialised)
        unInited.push_back(getName());

    for(auto component : components) {
        component->checkInitStatus(unInited);
    }
}

__host__ void PdfBase::recursiveSetNormalisation(fptype norm) const {
    host_normalisation[parameters] = norm;

    for(auto component : components) {
        component->recursiveSetNormalisation(norm);
    }
}

__host__ unsigned int PdfBase::registerParameter(Variable var) {
    static int unique_param = 0;

    if(find_in(parameterList, var))
        return static_cast<unsigned int>(var.getIndex());

    if(var.getIndex() < 0) {
        GOOFIT_DEBUG("{}: Registering p:{} for {}", getName(), unique_param, var.getName());
        var.setIndex(unique_param++);
    }

    parameterList.push_back(var);
    return static_cast<unsigned int>(var.getIndex());
}

__host__ void PdfBase::unregisterParameter(Variable var) {
    GOOFIT_DEBUG("{}: Removing {}", getName(), var.getName());

    for(PdfBase *comp : components) {
        comp->unregisterParameter(var);
    }

    var.setIndex(-1);
    // Once copies are used, this might able to be unregistred from a lower PDF only
    // For now, it gets completely cleared.
}

__host__ std::vector<Variable> PdfBase::getParameters() const {
    std::vector<Variable> ret;
    for(const Variable &param : parameterList)
        ret.push_back(param);

    for(const PdfBase *comp : components) {
        for(const Variable &sub_comp : comp->getParameters())
            if(!find_in(ret, sub_comp))
                ret.push_back(sub_comp);
    }

    return ret;
}

__host__ Variable *PdfBase::getParameterByName(std::string n) {
    for(Variable &p : parameterList) {
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
    for(const Observable &obs : observables)
        ret.push_back(obs);

    for(const PdfBase *comp : components) {
        for(const Observable &sub_comp : comp->getObservables())
            if(!find_in(ret, sub_comp))
                ret.push_back(sub_comp);
    }

    return ret;
}

__host__ unsigned int PdfBase::registerConstants(unsigned int amount) {
    if(totalConstants + amount >= maxParams)
        throw GooFit::GeneralError(
            "totalConstants {} + amount {} can not be more than {}", totalConstants, amount, maxParams);
    cIndex = totalConstants;
    totalConstants += amount;
    return cIndex;
}

void PdfBase::registerObservable(Observable obs) {
    if(find_in(observables, obs))
        return;

    GOOFIT_DEBUG("{}: Registering o:{} for {}", getName(), observables.size(), obs.getName());
    observables.push_back(obs);
}

__host__ void PdfBase::setIntegrationFineness(int i) {
    integrationBins = i;
    generateNormRange();
}

__host__ bool PdfBase::parametersChanged() const {
    return std::any_of(
        std::begin(parameterList), std::end(parameterList), [](const Variable &v) { return v.getChanged(); });
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
