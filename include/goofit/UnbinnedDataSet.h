#pragma once

#include "goofit/DataSet.h"

#include <initializer_list>
#include <map>
#include <vector>

#include <Eigen/Dense>

namespace GooFit {

class UnbinnedDataSet : public DataSet {
    // Class for unbinned datasets.

  public:
    using DataSet::addEvent;

    UnbinnedDataSet(Observable *var, std::string n = "");
    UnbinnedDataSet(std::vector<Observable *> &vars, std::string n = "");
    UnbinnedDataSet(std::set<Observable *> &vars, std::string n = "");
    UnbinnedDataSet(std::initializer_list<Observable *> vars, std::string n = "");

    ~UnbinnedDataSet() override = default;

    void addEvent() override;

    /// Get the value at a specific variable and event number
    fptype getValue(Observable *var, size_t idx) const;

    /// Set all the variables to the current event values
    void loadEvent(size_t idx);

    /// Set all entries to a constant value (note: this is kind of ugly)
    void setValueForAllEvents(Observable *var);

    /// Input an eigen matrix
    template <typename M>
    void from_matrix(const M &input, bool filter = false) {
        size_t optional_index = getNumEvents(); // Only used if index not included

        if(variables.size() != input.rows() && variables.size() != input.rows() + 1)
            throw GeneralError(
                "The wrong number of rows, expected {}, but matrix had {}", variables.size(), input.rows());

        for(int i = 0; i < input.cols(); i++) {     // Loop over events
            for(int j = 0; j < input.rows(); j++) { // Loop over variables
                variables.at(j)->setValue(input(j, i));
            }

            // Special override for counting variables (final param)
            if(variables.size() == input.rows() + 1)
                variables.at(input.rows())->setValue(optional_index++);

            if(!filter
               || std::all_of(std::begin(variables), std::end(variables), [](Observable *var) { return bool(*var); }))
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
