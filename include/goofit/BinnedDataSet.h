#pragma once

#include "goofit/DataSet.h"

#include <map>
#include <vector>

class BinnedDataSet : public DataSet {
    // Class for rectangularly binned datasets - every bin the same size.

public:
    using DataSet::addEvent;
    
    BinnedDataSet(Variable* var, std::string n = "");
    BinnedDataSet(std::vector<Variable*>& vars, std::string n = "");
    BinnedDataSet(std::set<Variable*>& vars, std::string n = "");
    BinnedDataSet(std::initializer_list<Variable*> vars, std::string n="");
    virtual ~BinnedDataSet() = default;

    virtual void addEvent() override;
    virtual void addWeightedEvent(double weight) override;

    fptype getBinContent(size_t bin) const {
        return binvalues[bin];
    }
    fptype getBinCenter(size_t ivar, size_t bin) const;
    fptype getBinCenter(Variable* var, size_t bin) const;
    
    size_t getBinNumber() const;
    fptype getBinVolume(size_t bin) const;
    fptype getBinError(size_t bin) const;
    
    size_t getNumBins() const;
    
    /// This includes weights
    fptype getNumWeightedEvents() const;

    void setBinContent(unsigned int bin, fptype value) {
        binvalues[bin] = value;
    }
    void setBinError(unsigned int bin, fptype error);

private:
    /// Lock all variable bin numbers
    void lockBins();
    
    std::vector<size_t> convertValuesToBins(const std::vector<fptype>& vals) const;
    size_t localToGlobal(const std::vector<size_t>& locals) const;
    std::vector<size_t> globalToLocal(size_t global) const;

    std::vector<fptype> binvalues;
    std::vector<fptype> binerrors;

};

