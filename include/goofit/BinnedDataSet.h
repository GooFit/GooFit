#pragma once

#include "goofit/DataSet.h"

#include <map>
#include <vector>

namespace GooFit {

class BinnedDataSet : public DataSet {
    // Class for rectangularly binned datasets - every bin the same size.

  public:
    using DataSet::addEvent;

    BinnedDataSet(Variable *var, std::string n = "");
    BinnedDataSet(std::vector<Variable *> &vars, std::string n = "");
    BinnedDataSet(std::set<Variable *> &vars, std::string n = "");
    BinnedDataSet(std::initializer_list<Variable *> vars, std::string n = "");
    ~BinnedDataSet() override = default;

    void addEvent() override;
    void addWeightedEvent(double weight) override;

    fptype getBinContent(size_t bin) const { return binvalues.at(bin); }
    fptype getBinCenter(size_t ivar, size_t bin) const;
    fptype getBinCenter(Variable *var, size_t bin) const;
    fptype getBinSize(size_t ivar) const;

    size_t getBinNumber() const;
    fptype getBinVolume(size_t bin) const;
    fptype getBinError(size_t bin) const;

    size_t getNumBins() const;

    /// This includes weights
    fptype getNumWeightedEvents() const;

    void setBinContent(unsigned int bin, fptype value) { binvalues.at(bin) = value; }
    void setBinError(unsigned int bin, fptype error);

  private:
    std::vector<size_t> convertValuesToBins(const std::vector<fptype> &vals) const;
    size_t localToGlobal(const std::vector<size_t> &locals) const;
    std::vector<size_t> globalToLocal(size_t global) const;

    /// Capture the number of bins on all variables
    void collectBins();

    std::vector<size_t> binsizes;
    std::vector<fptype> binvalues;
    std::vector<fptype> binerrors;
};

} // namespace GooFit
