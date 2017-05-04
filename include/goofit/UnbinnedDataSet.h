#pragma once

#include "goofit/DataSet.h"

#include <map>
#include <vector>
#include <initializer_list>

class UnbinnedDataSet : public DataSet {
    // Class for unbinned datasets.

public:
    using DataSet::addEvent;
    
    UnbinnedDataSet(Variable* var, std::string n = "");
    UnbinnedDataSet(std::vector<Variable*>& vars, std::string n = "");
    UnbinnedDataSet(std::set<Variable*>& vars, std::string n = "");
    UnbinnedDataSet(std::initializer_list<Variable*> vars, std::string n="");
    
    virtual ~UnbinnedDataSet();

    virtual void addEvent() override;
    int getNumEvents() const {
        return data.size();
    }
    fptype getValue(Variable* var, size_t idx) const;
    void loadEvent(size_t idx);
    void setValueForAllEvents(Variable* var);

private:
    std::vector<std::vector<fptype>> data;

};
