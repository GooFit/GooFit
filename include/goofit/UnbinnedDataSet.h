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

    UnbinnedDataSet(Variable *var, std::string n = "");
    UnbinnedDataSet(std::vector<Variable *> &vars, std::string n = "");
    UnbinnedDataSet(std::set<Variable *> &vars, std::string n = "");
    UnbinnedDataSet(std::initializer_list<Variable *> vars, std::string n = "");

    ~UnbinnedDataSet() override = default;

    void addEvent() override;

    /// Get the value at a specific variable and event number
    fptype getValue(Variable *var, size_t idx) const;

    /// Set all the variables to the current event values
    void loadEvent(size_t idx);

    /// Set all entries to a constant value (note: this is kind of ugly)
    void setValueForAllEvents(Variable *var);
    
    /// Input an eigen matrix
    template<typename M>
    void from_matrix(const M& input) {
        // Special override for counting variables (final param)
        if (data.size() == input.cols()+1) {
            size_t val = data[data.size()-1].size();
            for(int i=0; i<input.rows(); i++) {
                data[data.size()-1].push_back(val++);
            }
            
        } else if (data.size() != input.cols()) {
            throw GeneralError("The wrong number of rows, was {}, but matrix had {}", data.size(), input.cols());
        }
        for(int i=0; i<input.cols(); i++) {
            for(int j=0; j<input.rows(); j++) {
                data[i].push_back(input(j,i));
            }
        }
    }
    
    /// Produce an eigen Matrix
    template<typename M>
    M to_matrix() const {
        size_t rows = data.at(0).size();
        size_t columns = data.size();
        M mat(rows, columns);
        for (int i = 0; i < columns; i++) {
            for(int j=0; j<rows; j++)
                mat(j,i) = data[i].at(j);
        }
        return mat;
    }
    
  private:
    std::vector<std::vector<fptype>> data;
};
} // namespace GooFit
