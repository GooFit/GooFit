#pragma once

#include "goofit/DataSet.h"

#include <map>
#include <vector>
#include <initializer_list>

class UnbinnedDataSet : public DataSet {
    // Class for unbinned datasets.

public:
    using DataSet::DataSet;
    virtual ~UnbinnedDataSet();

    virtual void addEventVector(std::vector<fptype>& vals, fptype weight = 1);
    size_t getNumEvents() const {
        return data.size();
    }
    fptype getValue(Variable* var, int idx) const;
    void loadEvent(int idx);
    void setValueForAllEvents(Variable* var);

private:
    std::vector<std::map<Variable*, fptype>> data;

};
