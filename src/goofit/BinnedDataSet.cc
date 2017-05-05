#include "goofit/BinnedDataSet.h"
#include "goofit/Variable.h"
#include "goofit/Error.h"
#include "goofit/Log.h"

// Special constructor for one variable
BinnedDataSet::BinnedDataSet(Variable* var, std::string n)
    : DataSet(var, n) {
    lockBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::vector<Variable*>& vars, std::string n)
    : DataSet(vars, n) {
    lockBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::set<Variable*>& vars, std::string n)
    : DataSet(vars, n) {
    lockBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::initializer_list<Variable*> vars, std::string n)
: DataSet(vars, n) {
    lockBins();
    binvalues.resize(getNumBins());
}

void BinnedDataSet::addEvent() {
    size_t ibin = getBinNumber();
    binvalues.at(ibin) += 1;
    numEventsAdded++;
}

void BinnedDataSet::addWeightedEvent(double weight) {
    size_t ibin = getBinNumber();
    binvalues.at(ibin) += weight;
    numEventsAdded++;
    
}

void BinnedDataSet::lockBins() {
    for(Variable* var : variables)
        var->setLockedBins(true);
}

size_t BinnedDataSet::getBinNumber() const {
    std::vector<fptype> vals = getCurrentValues();
    std::vector<size_t> locals = convertValuesToBins(vals);
    return localToGlobal(locals);
}

size_t BinnedDataSet::localToGlobal(const std::vector<size_t>& locals) const {
    unsigned int priorMatrixSize = 1;
    unsigned int ret = 0;

    for(size_t i=0; i<variables.size(); i++) {
        unsigned int localBin = locals[i];
        ret += localBin * priorMatrixSize;
        priorMatrixSize *= variables[i]->getNumBins();
    }

    return ret;
}

std::vector<size_t> BinnedDataSet::globalToLocal(size_t global) const {
    std::vector<size_t> locals;

    // To convert global bin number to (x,y,z...) coordinates: For each dimension, take the mod
    // with the number of bins in that dimension. Then divide by the number of bins, in effect
    // collapsing so the grid has one fewer dimension. Rinse and repeat.
    
    for(size_t i=0; i<variables.size(); i++) {
        int localBin = global % variables[i]->getNumBins();
        locals.push_back(localBin);
        global /= variables[i]->getNumBins();
    }
    return locals;
}

fptype BinnedDataSet::getBinCenter(size_t ivar, size_t bin) const {
    std::vector<size_t> locals = globalToLocal(bin);
    size_t localBin = locals.at(ivar);
    
    fptype ret = variables[ivar]->getBinSize();
    ret       *= (localBin + 0.5);
    ret       += variables[ivar]->getLowerLimit();
    return ret;
}

fptype BinnedDataSet::getBinCenter(Variable* var, size_t bin) const {
    size_t ivar = indexOfVariable(var);
    return getBinCenter(ivar, bin);
}


fptype BinnedDataSet::getBinVolume(size_t bin) const {
    fptype ret = 1;

    for(size_t i=0; i<variables.size(); i++) {
        fptype step = variables[i]->getBinSize();
        ret *= step;
    }

    return ret;
}

fptype BinnedDataSet::getBinError(size_t bin) const {
    if(0 == binerrors.size())
        return sqrt(binvalues.at(bin));

    return binerrors.at(bin);
}

void BinnedDataSet::setBinError(unsigned int bin, fptype error) {
    if(0 == binerrors.size())
        binerrors.resize(binvalues.size());

    binerrors.at(bin) = error;
}

size_t BinnedDataSet::getNumBins() const {
    unsigned int ret = 1;

    for(size_t i=0; i<variables.size(); i++) {
        ret *= variables[i]->getNumBins();
    }

    return ret;
}

fptype BinnedDataSet::getNumEvents() const {
    fptype ret = 0;

    for(const fptype& bin : binvalues)
        ret += bin;

    return ret;
}

std::vector<size_t> BinnedDataSet::convertValuesToBins(const std::vector<fptype>& vals) const {
    if(vals.size() != variables.size())
        throw GooFit::GeneralError("Incorrect number of bins {} for {} variables", vals.size(), variables.size());

    std::vector<size_t> localBins;
    for(size_t i=0; i<variables.size(); i++) {
        fptype currval = vals[i];
        fptype betval = std::min(std::max(currval, variables[i]->getLowerLimit()),variables[i]->getUpperLimit());
        if(currval != betval)
            GOOFIT_INFO("Warning: Value {} outside {} range [{},{}] - clamping to {}",
                        currval, variables[i]->getName(), variables[i]->getLowerLimit(), variables[i]->getUpperLimit(), betval);
        fptype step = variables[i]->getBinSize();
        localBins.push_back( (size_t) floor((betval - variables[i]->getLowerLimit())/step));
    
    }
    

    return localBins;
}
