#pragma once

#include "goofit/DataSet.h"

#include <map>
#include <vector>

class UnbinnedDataSet : public DataSet {
    // Class for unbinned datasets.

public:
    UnbinnedDataSet(Variable* var, string n = "");
    UnbinnedDataSet(std::vector<Variable*>& vars, string n = "");
    UnbinnedDataSet(std::set<Variable*>& vars, string n = "");
    virtual ~UnbinnedDataSet();

    virtual void addEventVector(std::vector<fptype>& vals, fptype weight = 1);
    int getNumEvents() const {
        return data.size();
    }
    fptype getValue(Variable* var, int idx) const;
    void loadEvent(int idx);
    void setValueForAllEvents(Variable* var);

private:
    std::vector<std::map<Variable*, fptype>> data;

};
