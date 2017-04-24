#include "goofit/GlobalCudaDefines.h"
#include "goofit/PdfBase.h"
#include "goofit/Variable.h"
#include <algorithm>

fptype* dev_event_array;
fptype host_normalisation[maxParams];
fptype host_params[maxParams];
unsigned int host_indices[maxParams];
int host_callnumber = 0;
int totalParams = 0;
int totalConstants = 1; // First constant is reserved for number of events.
std::map<Variable*, std::set<PdfBase*>> variableRegistry;

PdfBase::PdfBase(Variable* x, std::string n)
    : name(n) { // Special-case PDFs should set to false.
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
    host_normalisation[parameters] = norm;

    for(unsigned int i = 0; i < components.size(); ++i) {
        components[i]->recursiveSetNormalisation(norm);
    }
}

__host__ unsigned int PdfBase::registerParameter(Variable* var) {
    if(!var) {
        std::cout << "Error: Attempt to register null Variable with "
                  << getName()
                  << ", aborting.\n";
        assert(var);
        exit(1);
    }

    if(std::find(parameterList.begin(), parameterList.end(), var) != parameterList.end())
        return (unsigned int) var->getIndex();

    parameterList.push_back(var);
    variableRegistry[var].insert(this);

    if(0 > var->getIndex()) {
        unsigned int unusedIndex = 0;

        while(true) {
            bool canUse = true;

            for(std::map<Variable*, std::set<PdfBase*>>::iterator p = variableRegistry.begin(); p != variableRegistry.end(); ++p) {
                if(unusedIndex != (*p).first->index)
                    continue;

                canUse = false;
                break;
            }

            if(canUse)
                break;

            unusedIndex++;
        }

        var->index = unusedIndex;
    }

    return (unsigned int) var->getIndex();
}

__host__ void PdfBase::unregisterParameter(Variable* var) {
    if(!var)
        return;

    parIter pos = std::find(parameterList.begin(), parameterList.end(), var);

    if(pos != parameterList.end())
        parameterList.erase(pos);

    variableRegistry[var].erase(this);

    if(0 == variableRegistry[var].size())
        var->index = -1;

    for(unsigned int i = 0; i < components.size(); ++i) {
        components[i]->unregisterParameter(var);
    }
}

__host__ void PdfBase::getParameters(parCont& ret) const {
    for(parConstIter p = parameterList.begin(); p != parameterList.end(); ++p) {
        if(std::find(ret.begin(), ret.end(), (*p)) != ret.end())
            continue;

        ret.push_back(*p);
    }

    for(unsigned int i = 0; i < components.size(); ++i) {
        components[i]->getParameters(ret);
    }
}

__host__ Variable* PdfBase::getParameterByName(std::string n) const {
    for(parConstIter p = parameterList.begin(); p != parameterList.end(); ++p) {
        if((*p)->name == n)
            return (*p);
    }

    for(unsigned int i = 0; i < components.size(); ++i) {
        Variable* cand = components[i]->getParameterByName(n);

        if(cand)
            return cand;
    }

    return 0;
}

__host__ void PdfBase::getObservables(std::vector<Variable*>& ret) const {
    for(obsConstIter p = obsCBegin(); p != obsCEnd(); ++p) {
        if(std::find(ret.begin(), ret.end(), *p) != ret.end())
            continue;

        ret.push_back(*p);
    }

    for(unsigned int i = 0; i < components.size(); ++i) {
        components[i]->getObservables(ret);
    }
}

__host__ unsigned int PdfBase::registerConstants(unsigned int amount) {
    assert(totalConstants + amount < maxParams);
    cIndex = totalConstants;
    totalConstants += amount;
    return cIndex;
}

void PdfBase::registerObservable(Variable* obs) {
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
    for(Variable* v : parameterList) {
        if(v->changed())
            return true;
    }
    return false;
}

__host__ void PdfBase::storeParameters() const {
    parCont params;
    getParameters(params);

    if(!cachedParams)
        cachedParams = new fptype[params.size()];

    int counter = 0;

    for(Variable* v : params) {
        cachedParams[counter++] = host_params[v->index];
    }
}

__host__ void PdfBase::setNumPerTask(PdfBase* p, const int& c) {
    if(!p)
        return;

    m_iEventsPerTask = c;
}

