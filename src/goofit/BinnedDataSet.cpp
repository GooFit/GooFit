#include <goofit/BinnedDataSet.h>
#include <goofit/Error.h>
#include <goofit/Log.h>
#include <goofit/Variable.h>

#include <functional>
#include <numeric>

namespace GooFit {

// Special constructor for one variable
BinnedDataSet::BinnedDataSet(const Observable &var, std::string n)
    : DataSet(var, n) {
    collectBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(const std::vector<Observable> &vars, std::string n)
    : DataSet(vars, n) {
    collectBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(const std::set<Observable> &vars, std::string n)
    : DataSet(vars, n) {
    collectBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::initializer_list<Observable> vars, std::string n)
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

    for(const Observable &var : observables)
        binsizes.push_back(var.getNumBins());
}

size_t BinnedDataSet::getBinNumber() const {
    std::vector<fptype> vals   = getCurrentValues();
    std::vector<size_t> locals = convertValuesToBins(vals);
    return localToGlobal(locals);
}

size_t BinnedDataSet::localToGlobal(const std::vector<size_t> &locals) const {
    unsigned int priorMatrixSize = 1;
    unsigned int ret             = 0;

    for(size_t i = 0; i < observables.size(); i++) {
        size_t localBin = locals[i];
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

    for(size_t i = 0; i < observables.size(); i++) {
        size_t localBin = global % binsizes[i];
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
    ret += observables[ivar].getLowerLimit();
    return ret;
}

fptype BinnedDataSet::getBinCenter(const Observable &var, size_t bin) const {
    size_t ivar = indexOfVariable(var);
    return getBinCenter(ivar, bin);
}

fptype BinnedDataSet::getBinSize(size_t ivar) const {
    return (observables.at(ivar).getUpperLimit() - observables[ivar].getLowerLimit()) / binsizes[ivar];
}

fptype BinnedDataSet::getBinVolume(size_t bin) const {
    fptype ret = 1;

    for(size_t i = 0; i < observables.size(); i++) {
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
    if(vals.size() != observables.size())
        throw GooFit::GeneralError("Incorrect number of bins {} for {} variables", vals.size(), observables.size());

    std::vector<size_t> localBins;
    for(size_t i = 0; i < observables.size(); i++) {
        fptype currval = vals[i];
        fptype betval  = std::min(std::max(currval, observables[i].getLowerLimit()), observables[i].getUpperLimit());
        if(currval != betval)
            GOOFIT_INFO("Warning: Value {} outside {} range [{},{}] - clamping to {}",
                        currval,
                        observables[i].getName(),
                        observables[i].getLowerLimit(),
                        observables[i].getUpperLimit(),
                        betval);
        localBins.push_back(static_cast<size_t>(floor((betval - observables[i].getLowerLimit()) / getBinSize(i))));
    }

    return localBins;
}
} // namespace GooFit
