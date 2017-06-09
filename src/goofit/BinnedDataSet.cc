#include "goofit/BinnedDataSet.h"
#include "goofit/Error.h"
#include "goofit/Log.h"
#include "goofit/Variable.h"

#include <functional>
#include <numeric>

namespace GooFit {

// Special constructor for one variable
BinnedDataSet::BinnedDataSet(Variable *var, std::string n)
    : DataSet(var, n) {
    collectBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::vector<Variable *> &vars, std::string n)
    : DataSet(vars, n) {
    collectBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::set<Variable *> &vars, std::string n)
    : DataSet(vars, n) {
    collectBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::initializer_list<Variable *> vars, std::string n)
    : DataSet(vars, n) {
    collectBins();
    binvalues.resize(getNumBins());
}

void BinnedDataSet::addEvent() {
    checkAllVars();
    size_t ibin = getBinNumber();
    binvalues.at(ibin) += 1;
    numEventsAdded++;
}

void BinnedDataSet::addWeightedEvent(double weight) {
    checkAllVars();
    size_t ibin = getBinNumber();
    binvalues.at(ibin) += weight;
    numEventsAdded++;
}

void BinnedDataSet::collectBins() {
    // Not really intended to be run multiple times, but just in case
    binsizes.clear();

    for(Variable *var : variables)
        binsizes.push_back(var->getNumBins());
}

size_t BinnedDataSet::getBinNumber() const {
    std::vector<fptype> vals   = getCurrentValues();
    std::vector<size_t> locals = convertValuesToBins(vals);
    return localToGlobal(locals);
}

size_t BinnedDataSet::localToGlobal(const std::vector<size_t> &locals) const {
    unsigned int priorMatrixSize = 1;
    unsigned int ret             = 0;

    for(size_t i = 0; i < variables.size(); i++) {
        unsigned int localBin = locals[i];
        ret += localBin * priorMatrixSize;
        priorMatrixSize *= binsizes[i];
    }

    return ret;
}

std::vector<size_t> BinnedDataSet::globalToLocal(size_t global) const {
    std::vector<size_t> locals;

    // To convert global bin number to (x,y,z...) coordinates: For each dimension, take the mod
    // with the number of bins in that dimension. Then divide by the number of bins, in effect
    // collapsing so the grid has one fewer dimension. Rinse and repeat.

    for(size_t i = 0; i < variables.size(); i++) {
        int localBin = global % binsizes[i];
        locals.push_back(localBin);
        global /= binsizes[i];
    }
    return locals;
}

fptype BinnedDataSet::getBinCenter(size_t ivar, size_t bin) const {
    std::vector<size_t> locals = globalToLocal(bin);
    size_t localBin            = locals.at(ivar);

    fptype ret = getBinSize(ivar);
    ret *= (localBin + 0.5);
    ret += variables[ivar]->getLowerLimit();
    return ret;
}

fptype BinnedDataSet::getBinCenter(Variable *var, size_t bin) const {
    size_t ivar = indexOfVariable(var);
    return getBinCenter(ivar, bin);
}

fptype BinnedDataSet::getBinSize(size_t ivar) const {
    return (variables.at(ivar)->getUpperLimit() - variables[ivar]->getLowerLimit()) / binsizes[ivar];
}

fptype BinnedDataSet::getBinVolume(size_t bin) const {
    fptype ret = 1;

    for(size_t i = 0; i < variables.size(); i++) {
        ret *= getBinSize(i);
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
    return std::accumulate(std::begin(binsizes), std::end(binsizes), 1, std::multiplies<size_t>());
}

fptype BinnedDataSet::getNumWeightedEvents() const {
    return std::accumulate(std::begin(binvalues), std::end(binvalues), 0);
}

std::vector<size_t> BinnedDataSet::convertValuesToBins(const std::vector<fptype> &vals) const {
    if(vals.size() != variables.size())
        throw GooFit::GeneralError("Incorrect number of bins {} for {} variables", vals.size(), variables.size());

    std::vector<size_t> localBins;
    for(size_t i = 0; i < variables.size(); i++) {
        fptype currval = vals[i];
        fptype betval  = std::min(std::max(currval, variables[i]->getLowerLimit()), variables[i]->getUpperLimit());
        if(currval != betval)
            GOOFIT_INFO("Warning: Value {} outside {} range [{},{}] - clamping to {}",
                        currval,
                        variables[i]->getName(),
                        variables[i]->getLowerLimit(),
                        variables[i]->getUpperLimit(),
                        betval);
        localBins.push_back(static_cast<size_t>(floor((betval - variables[i]->getLowerLimit()) / getBinSize(i))));
    }

    return localBins;
}
} // namespace GooFit
