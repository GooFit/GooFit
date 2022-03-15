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
    auto getBinContent(size_t bin) const -> fptype { return binvalues.at(bin); }
    auto getBinCenter(size_t ivar, size_t bin) const -> fptype;
    auto getBinCenter(const Observable &var, size_t bin) const -> fptype;

    /// Get the size of a bin
    auto getBinSize(size_t ivar) const -> fptype;

    auto getBinNumber() const -> size_t;
    auto getBinVolume(size_t bin) const -> fptype;
    auto getBinError(size_t bin) const -> fptype;

    /// Get the size of each diminsion
    auto getDiminsions() const -> std::vector<size_t> { return binsizes; }

    auto getNumBins() const -> size_t;

    /// This includes weights
    auto getNumWeightedEvents() const -> fptype;

    void setBinContent(unsigned int bin, fptype value) { binvalues.at(bin) = value; }
    void setBinError(unsigned int bin, fptype value);

  private:
    auto convertValuesToBins(const std::vector<fptype> &vals) const -> std::vector<size_t>;
    auto localToGlobal(const std::vector<size_t> &locals) const -> size_t;
    auto globalToLocal(size_t global) const -> std::vector<size_t>;

    /// Capture the number of bins on all variables
    void collectBins();

    std::vector<size_t> binsizes;
    std::vector<fptype> binvalues;
    std::vector<fptype> binerrors;
};

} // namespace GooFit
