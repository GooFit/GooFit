#include "goofit/GlobalCudaDefines.h"
#include "goofit/PdfBase.h"
#include "goofit/Variable.h"
#include "goofit/Color.h"
#include "goofit/Error.h"
#include "goofit/Log.h"

#include <algorithm>
#include <set>

#include "goofit/FitManager.h"
#include "goofit/BinnedDataSet.h"
#include "goofit/UnbinnedDataSet.h"

namespace GooFit {


fptype* dev_event_array;
fptype host_normalisations[maxParams];
fptype host_parameters[maxParams];
unsigned int host_constants[maxParams];
fptype host_observables[maxParams];

int host_callnumber = 0;

int totalParameters = 0;
int totalConstants = 0; // First constant is reserved for number of events.
int totalObservables = 0;
int totalNormalisations = 0;

std::map<Variable*, std::set<PdfBase*>> variableRegistry;

PdfBase::PdfBase(Variable* x, std::string n)
    : name(n), parametersIdx(0) { // Special-case PDFs should set to false.
    if(x)
        registerObservable(x);
}

__host__ void PdfBase::checkInitStatus(std::vector<std::string>& unInited) const {
    if(!properlyInitialised)
        unInited.push_back(getName());

    for(unsigned int i = 0; i < components.size(); ++i) {
        components[i]->checkInitStatus(unInited);
    }
}

__host__ void PdfBase::recursiveSetNormalisation(fptype norm) const {
    host_normalisations[normalIdx + 1] = norm;

    for(unsigned int i = 0; i < components.size(); ++i) {
        components[i]->recursiveSetNormalisation(norm);
    }
}

__host__ unsigned int PdfBase::registerParameter(Variable* var) {
    if(var == nullptr)
        throw GooFit::GeneralError("{}: Can not register a nullptr", getName());

    //we need to always add the variable to our internal parameter list
    parametersList.push_back(var);

    if(std::find(parametersList.begin(), parametersList.end(), var) != parametersList.end())
        return (unsigned int) var->getIndex();

    //parametersList.push_back(var);
    variableRegistry[var].insert(this);

    if(0 > var->getIndex()) {
        unsigned int unusedIndex = 0;

        while(true) {
            bool canUse = true;

            for(std::map<Variable*, std::set<PdfBase*>>::iterator p = variableRegistry.begin(); p != variableRegistry.end(); ++p) {
                if(unusedIndex != (*p).first->getIndex())
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

    return (unsigned int) var->getIndex();
}

__host__ void PdfBase::unregisterParameter(Variable* var) {
    if(var == nullptr)
        return;

    GOOFIT_DEBUG("{}: Removing {}", getName(), var->getName());
    
    auto pos = std::find(parametersList.begin(), parametersList.end(), var);

    if(pos != parametersList.end())
        parametersList.erase(pos);

    variableRegistry[var].erase(this);

    if(0 == variableRegistry[var].size())
        var->setIndex(-1);

    for(PdfBase* comp : components) {
        comp->unregisterParameter(var);
    }
}


__host__ std::vector<Variable*> PdfBase::getParameters() const {
    
    std::vector<Variable*> ret = parametersList;
    
    for(const PdfBase* comp : components) {
        for(Variable* sub_comp : comp->getParameters())
            if(std::find(std::begin(ret), std::end(ret), sub_comp)==std::end(ret))
                ret.push_back(sub_comp);
    }
    
    return ret;
}


__host__ Variable* PdfBase::getParameterByName(std::string n) const {
    for(Variable* p : parametersList) {
        if(p->getName() == n)
            return p;
    }

    for(unsigned int i = 0; i < components.size(); ++i) {
        Variable* cand = components[i]->getParameterByName(n);

        if(cand)
            return cand;
    }

    return 0;
}

__host__ std::vector<Variable*> PdfBase::getObservables() const {
    std::vector<Variable*> ret = observablesList;
    
    for(const PdfBase* comp : components) {
        for(Variable* sub_comp : comp->getObservables())
            if(std::find(std::begin(ret), std::end(ret), sub_comp)==std::end(ret))
                ret.push_back(sub_comp);
    }
    
    return ret;
}

__host__ unsigned int PdfBase::registerConstants(unsigned int amount) {
    if(totalConstants + amount >= maxParams)
        throw GooFit::GeneralError("totalConstants {} + amount {} can not be more than {}", totalConstants, amount, maxParams);
    cIndex = totalConstants;
    totalConstants += amount;
    return cIndex;
}

void PdfBase::registerObservable(Variable* obs) {
    if(!obs)
        return;

    if(find(observablesList.begin(), observablesList.end(), obs) != observablesList.end())
        return;

    observablesList.push_back(obs);
}

__host__ void PdfBase::setIntegrationFineness(int i) {
    integrationBins = i;
    generateNormRange();
}

__host__ bool PdfBase::parametersChanged() const {
    return std::any_of(std::begin(parametersList), std::end(parametersList), [](Variable* v){return v->getChanged();});
}

__host__ void PdfBase::setNumPerTask(PdfBase* p, const int& c) {
    if(!p)
        return;

    m_iEventsPerTask = c;
}

__host__ ROOT::Minuit2::FunctionMinimum PdfBase::fitTo(DataSet *data) {
    setData(data);
    FitManager fitter{this};
    return fitter.fit();
}
} // namespace GooFit

