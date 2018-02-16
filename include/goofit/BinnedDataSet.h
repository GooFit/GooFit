#pragma once

#include <map>
#include <vector>

#include <goofit/DataSet.h>

namespace GooFit {

class BinnedDataSet : public DataSet {
    // Class for rectangularly binned datasets - every bin the same size.

  public:
    using DataSet::addEvent;

    BinnedDataSet(const Observable &var, std::string n = "");
    BinnedDataSet(const std::vector<Observable> &vars, std::string n = "");
    BinnedDataSet(const std::set<Observable> &vars, std::string n = "");
    BinnedDataSet(std::initializer_list<Observable> vars, std::string n = "");
    ~BinnedDataSet() override = default;

    void addEvent() override;
    void addWeightedEvent(double weight) override;

    /// Get the content of a bin
    fptype getBinContent(size_t bin) const { return binvalues.at(bin); }
    fptype getBinCenter(size_t ivar, size_t bin) const;
    fptype getBinCenter(const Observable &var, size_t bin) const;

    /// Get the size of a bin
    fptype getBinSize(size_t ivar) const;

    size_t getBinNumber() const;
    fptype getBinVolume(size_t bin) const;
    fptype getBinError(size_t bin) const;

    /// Get the size of each diminsion
    std::vector<size_t> getDiminsions() const { return binsizes; }

    size_t getNumBins() const;

    /// This includes weights
    fptype getNumWeightedEvents() const;

    void setBinContent(unsigned int bin, fptype value) { binvalues.at(bin) = value; }
    void setBinError(unsigned int bin, fptype value);

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
