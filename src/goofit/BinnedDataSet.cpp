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

auto BinnedDataSet::getBinNumber() const -> size_t {
    std::vector<fptype> vals   = getCurrentValues();
    std::vector<size_t> locals = convertValuesToBins(vals);
    return localToGlobal(locals);
}

auto BinnedDataSet::localToGlobal(const std::vector<size_t> &locals) const -> size_t {
    unsigned int priorMatrixSize = 1;
    unsigned int ret             = 0;

    for(size_t i = 0; i < observables.size(); i++) {
        size_t localBin = locals[i];
        ret += localBin * priorMatrixSize;
        priorMatrixSize *= binsizes[i];
    }

    return ret;
}

auto BinnedDataSet::globalToLocal(size_t global) const -> std::vector<size_t> {
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

auto BinnedDataSet::getBinCenter(size_t ivar, size_t bin) const -> fptype {
    std::vector<size_t> locals = globalToLocal(bin);
    size_t localBin            = locals.at(ivar);

    fptype ret = getBinSize(ivar);
    ret *= (localBin + 0.5);
    ret += observables[ivar].getLowerLimit();
    return ret;
}

auto BinnedDataSet::getBinCenter(const Observable &var, size_t bin) const -> fptype {
    size_t ivar = indexOfVariable(var);
    return getBinCenter(ivar, bin);
}

auto BinnedDataSet::getBinSize(size_t ivar) const -> fptype {
    return (observables.at(ivar).getUpperLimit() - observables[ivar].getLowerLimit()) / binsizes[ivar];
}

auto BinnedDataSet::getBinVolume(size_t bin) const -> fptype {
    fptype ret = 1;

    for(size_t i = 0; i < observables.size(); i++) {
        ret *= getBinSize(i);
    }

    return ret;
}

auto BinnedDataSet::getBinError(size_t bin) const -> fptype {
    if(0 == binerrors.size())
        return sqrt(binvalues.at(bin));

    return binerrors.at(bin);
}

void BinnedDataSet::setBinError(unsigned int bin, fptype error) {
    if(0 == binerrors.size())
        binerrors.resize(binvalues.size());

    binerrors.at(bin) = error;
}

auto BinnedDataSet::getNumBins() const -> size_t {
    return std::accumulate(std::begin(binsizes), std::end(binsizes), 1, std::multiplies<size_t>());
}

auto BinnedDataSet::getNumWeightedEvents() const -> fptype {
    return std::accumulate(std::begin(binvalues), std::end(binvalues), 0);
}

auto BinnedDataSet::convertValuesToBins(const std::vector<fptype> &vals) const -> std::vector<size_t> {
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
