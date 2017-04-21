#include "goofit/BinnedDataSet.h"
#include "goofit/Variable.h"
#include <cassert>

// Special constructor for one variable
BinnedDataSet::BinnedDataSet(Variable* var, string n)
    : DataSet(var, n) {
    cacheNumBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::vector<Variable*>& vars, string n)
    : DataSet(vars, n) {
    cacheNumBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::BinnedDataSet(std::set<Variable*>& vars, string n)
    : DataSet(vars, n) {
    cacheNumBins();
    binvalues.resize(getNumBins());
}

BinnedDataSet::~BinnedDataSet() {}

void BinnedDataSet::addEventVector(std::vector<fptype>& vals, fptype weight) {
    numEventsAdded++;
    std::vector<unsigned int> localBins = convertValuesToBins(vals);
    unsigned int bin = localToGlobal(localBins);

    if(bin >= binvalues.size()) {
        std::cout << "Bad bin number "
                  << bin << " / "
                  << binvalues.size()
                  << " from input values ";

        for(unsigned int i = 0; i < vals.size(); ++i) {
            std::cout << vals[i] << " ";
        }

        std::cout << "\n";
        assert(bin < binvalues.size());
    }

    binvalues[bin] += weight;
}

void BinnedDataSet::cacheNumBins() {
    for(varConstIt v = varsBegin(); v != varsEnd(); ++v) {
        cachedNumBins[*v] = (*v)->numbins;
    }
}

unsigned int BinnedDataSet::getBinNumber() const {
    vector<fptype> vals = getCurrentValues();
    vector<unsigned int> locals = convertValuesToBins(vals);
    return localToGlobal(locals);
}

unsigned int BinnedDataSet::localToGlobal(std::vector<unsigned int>& locals) const {
    unsigned int priorMatrixSize = 1;
    unsigned int ret = 0;

    for(varConstIt v = varsBegin(); v != varsEnd(); ++v) {
        unsigned int localBin = locals[indexOfVariable(*v)];
        ret += localBin * priorMatrixSize;
        priorMatrixSize *= cachedNumBins.at(*v); // Use 'at' to preserve const-ness.
    }

    return ret;
}

void BinnedDataSet::globalToLocal(std::vector<unsigned int>& locals, unsigned int global) const {
    locals.clear();

    // To convert global bin number to (x,y,z...) coordinates: For each dimension, take the mod
    // with the number of bins in that dimension. Then divide by the number of bins, in effect
    // collapsing so the grid has one fewer dimension. Rinse and repeat.
    for(varConstIt v = varsBegin(); v != varsEnd(); ++v) {
        int localBin = global % cachedNumBins.at(*v);
        locals.push_back(localBin);
        global /= cachedNumBins.at(*v);
    }
}


fptype BinnedDataSet::getBinCenter(Variable* var, unsigned int bin) const {
    std::vector<unsigned int> locals;
    globalToLocal(locals, bin);
    unsigned int varIndex = indexOfVariable(var);
    unsigned int localBin = locals[varIndex];
    fptype ret = var->upperlimit;
    ret       -= var->lowerlimit;
    ret       /= cachedNumBins.at(var);
    ret       *= (localBin + 0.5);
    ret       += var->lowerlimit;
    return ret;
}

fptype BinnedDataSet::getBinVolume(unsigned int bin) const {
    fptype ret = 1;

    for(varConstIt v = varsBegin(); v != varsEnd(); ++v) {
        fptype step = (*v)->upperlimit;
        step       -= (*v)->lowerlimit;
        step       /= cachedNumBins.at(*v);
        ret *= step;
    }

    return ret;
}

fptype BinnedDataSet::getBinError(unsigned int bin) const {
    if(0 == binerrors.size())
        return SQRT(binvalues[bin]);

    assert(bin < binerrors.size());
    return binerrors[bin];
}

void BinnedDataSet::setBinError(unsigned int bin, fptype error) {
    if(0 == binerrors.size())
        binerrors.resize(binvalues.size());

    assert(bin < binerrors.size());
    binerrors[bin] = error;
}

unsigned int BinnedDataSet::getNumBins() const {
    unsigned int ret = 1;

    for(varConstIt v = varsBegin(); v != varsEnd(); ++v) {
        ret *= cachedNumBins.at(*v);
    }

    return ret;
}

fptype BinnedDataSet::getNumEvents() const {
    fptype ret = 0;

    for(std::vector<fptype>::const_iterator bin = binvalues.begin(); bin != binvalues.end(); ++bin)
        ret += (*bin);

    return ret;
}

vector<unsigned int> BinnedDataSet::convertValuesToBins(const vector<fptype>& vals) const {
    std::vector<unsigned int> localBins;
    varConstIt currVar = varsBegin();

    for(unsigned int i = 0; i < vals.size(); ++i) {
        assert(currVar != varsEnd());
        double currVal = vals[i];

        if(currVal < (*currVar)->lowerlimit) {
            std::cout << "Warning: Value "
                      << currVal
                      << " less than minimum "
                      << (*currVar)->lowerlimit
                      << " for "
                      << (*currVar)->name
                      << "; clamping to minimum.\n";
            currVal = (*currVar)->lowerlimit;
        }

        if(currVal > (*currVar)->upperlimit) {
            std::cout << "Warning: Value "
                      << currVal
                      << " more than maximum "
                      << (*currVar)->upperlimit
                      << " for "
                      << (*currVar)->name
                      << "; clamping to maximum.\n";
            currVal = (*currVar)->upperlimit;
        }

        fptype step = (*currVar)->upperlimit;
        step       -= (*currVar)->lowerlimit;
        step       /= cachedNumBins.at(*currVar);

        fptype curr = currVal;
        curr -= (*currVar)->lowerlimit;
        curr /= step;
        localBins.push_back((unsigned int) floor(curr));

        ++currVar;
    }

    return localBins;
}
