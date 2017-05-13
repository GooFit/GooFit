#pragma once

#include "goofit/DataSet.h"

#include <map>
#include <vector>

class BinnedDataSet : public DataSet {
    // Class for rectangularly binned datasets - every bin the same size.

public:
    BinnedDataSet(Variable* var, std::string n = "");
    BinnedDataSet(std::vector<Variable*>& vars, std::string n = "");
    BinnedDataSet(std::set<Variable*>& vars, std::string n = "");
    BinnedDataSet(std::initializer_list<Variable*> vars, std::string n="");
    virtual ~BinnedDataSet();

    virtual void addEventVector(std::vector<fptype>& vals, fptype weight = 1);

    fptype getBinContent(unsigned int bin) const {
        return binvalues[bin];
    }
    fptype getBinCenter(Variable* var, unsigned int bin) const;
    unsigned int getBinNumber() const;
    fptype getBinVolume(unsigned int bin) const;
    fptype getBinError(unsigned int bin) const;
    unsigned int getNumBins() const;
    size_t getNumEvents() const;

    void setBinContent(unsigned int bin, fptype value) {
        binvalues[bin] = value;
    }
    void setBinError(unsigned int bin, fptype error);

private:
    void cacheNumBins();
    std::vector<unsigned int> convertValuesToBins(const std::vector<fptype>& vals) const;
    unsigned int localToGlobal(std::vector<unsigned int>& locals) const;
    void globalToLocal(std::vector<unsigned int>& locals, unsigned int global) const;

    std::vector<fptype> binvalues;
    std::vector<fptype> binerrors;
    
    /// Store these numbers in case they change on the user end - vast confusion possible.
    std::map<Variable*, int> cachedNumBins;
};

