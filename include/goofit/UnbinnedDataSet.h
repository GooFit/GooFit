#pragma once

#include <goofit/DataSet.h>

#include <initializer_list>
#include <map>
#include <vector>

#include <Eigen/Dense>

namespace GooFit {

class UnbinnedDataSet : public DataSet {
    // Class for unbinned datasets.

  public:
    using DataSet::addEvent;

    UnbinnedDataSet(const Observable &var, std::string n = "");
    UnbinnedDataSet(const std::vector<Observable> &vars, std::string n = "");
    UnbinnedDataSet(const std::set<Observable> &vars, std::string n = "");
    UnbinnedDataSet(std::initializer_list<Observable> vars, std::string n = "");

    ~UnbinnedDataSet() override = default;

    void addEvent() override;

    /// Replace the current dataset with a grid
    void fillWithGrid();

    /// Get the value at a specific variable and event number
    fptype getValue(const Observable &var, size_t idx) const;

    /// Set all the variables to the current event values
    void loadEvent(size_t idx);

    /// Set all entries to a constant value (note: this is kind of ugly)
    void setValueForAllEvents(const Observable &var);

    /// Input an eigen matrix
    template <typename M>
    void from_matrix(const M &input, bool filter = false) {
        size_t optional_index = getNumEvents(); // Only used if index not included

        if(observables.size() != input.rows() && observables.size() != input.rows() + 1)
            throw GeneralError(
                "The wrong number of rows, expected {}, but matrix had {}", observables.size(), input.rows());

        for(int i = 0; i < input.cols(); i++) {     // Loop over events
            for(int j = 0; j < input.rows(); j++) { // Loop over variables
                observables.at(j).setValue(input(j, i));
            }

            // Special override for counting variables (final param)
            if(observables.size() == input.rows() + 1)
                observables.at(input.rows()).setValue(optional_index++);

            if(!filter
               || std::all_of(std::begin(observables), std::end(observables), [](Observable var) { return bool(var); }))
                addEvent();
        }
    }

    /// Produce an eigen Matrix
    template <typename M>
    M to_matrix() const {
        size_t rows = data.size();
        size_t cols = data.at(0).size();
        M mat(rows, cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++)
                mat(i, j) = data[i].at(j);
        }
        return mat;
    }

  private:
    std::vector<std::vector<fptype>> data;
};
} // namespace GooFit
