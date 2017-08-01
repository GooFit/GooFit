#include "goofit/GlobalCudaDefines.h"
#include "goofit/Color.h"
#include "goofit/Error.h"
#include "goofit/Log.h"
#include "goofit/PdfBase.h"
#include "goofit/Variable.h"

#include <algorithm>
#include <set>
#include <utility>

#include "goofit/BinnedDataSet.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"

namespace GooFit {

fptype *dev_event_array;
fptype host_normalisation[maxParams];
fptype host_params[maxParams];
unsigned int host_indices[maxParams];

int host_callnumber = 0;
int totalParams     = 0;
int totalConstants  = 1; // First constant is reserved for number of events.
std::map<Variable *, std::set<PdfBase *>> variableRegistry;

PdfBase::PdfBase(Variable *x, std::string n)
    : name(std::move(n)) { // Special-case PDFs should set to false.
    if(x)
        registerObservable(x);
}

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

__host__ unsigned int PdfBase::registerParameter(Variable *var) {
    if(var == nullptr)
        throw GooFit::GeneralError("{}: Can not register a nullptr", getName());

    if(std::find(parameterList.begin(), parameterList.end(), var) != parameterList.end())
        return static_cast<unsigned int>(var->getIndex());

    parameterList.push_back(var);
    variableRegistry[var].insert(this);

    if(0 > var->getIndex()) {
        unsigned int unusedIndex = 0;

        while(true) {
            bool canUse = true;

            for(auto &p : variableRegistry) {
                if(unusedIndex != p.first->getIndex())
                    continue;

                canUse = false;
                break;
            }

            if(canUse)
                break;

            unusedIndex++;
        }

        GOOFIT_DEBUG("{}: Registering {} for {}", getName(), unusedIndex, var->getName());
        var->setIndex(unusedIndex);
    }

    return static_cast<unsigned int>(var->getIndex());
}

__host__ void PdfBase::unregisterParameter(Variable *var) {
    if(var == nullptr)
        return;

    GOOFIT_DEBUG("{}: Removing {}", getName(), var->getName());

    auto pos = std::find(parameterList.begin(), parameterList.end(), var);

    if(pos != parameterList.end())
        parameterList.erase(pos);

    variableRegistry[var].erase(this);

    if(0 == variableRegistry[var].size())
        var->setIndex(-1);

    for(PdfBase *comp : components) {
        comp->unregisterParameter(var);
    }
}

__host__ std::vector<Variable *> PdfBase::getParameters() const {
    std::vector<Variable *> ret = parameterList;

    for(const PdfBase *comp : components) {
        for(Variable *sub_comp : comp->getParameters())
            if(std::find(std::begin(ret), std::end(ret), sub_comp) == std::end(ret))
                ret.push_back(sub_comp);
    }

    return ret;
}

__host__ Variable *PdfBase::getParameterByName(std::string n) const {
    for(Variable *p : parameterList) {
        if(p->getName() == n)
            return p;
    }

    for(auto component : components) {
        Variable *cand = component->getParameterByName(n);

        if(cand)
            return cand;
    }

    return nullptr;
}

__host__ std::vector<Variable *> PdfBase::getObservables() const {
    std::vector<Variable *> ret = observables;

    for(const PdfBase *comp : components) {
        for(Variable *sub_comp : comp->getObservables())
            if(std::find(std::begin(ret), std::end(ret), sub_comp) == std::end(ret))
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

void PdfBase::registerObservable(Variable *obs) {
    if(!obs)
        return;

    if(find(observables.begin(), observables.end(), obs) != observables.end())
        return;

    observables.push_back(obs);
}

__host__ void PdfBase::setIntegrationFineness(int i) {
    integrationBins = i;
    generateNormRange();
}

__host__ bool PdfBase::parametersChanged() const {
    return std::any_of(std::begin(parameterList), std::end(parameterList), [](Variable *v) { return v->getChanged(); });
}

__host__ void PdfBase::setNumPerTask(PdfBase *p, const int &c) {
    if(!p)
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
