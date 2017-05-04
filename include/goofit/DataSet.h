#pragma once

#include "goofit/GlobalCudaDefines.h"

#include <vector>
#include <set>
#include <string>
#include <initializer_list>

#include "goofit/Variable.h"

class DataSet {
public:
    DataSet(Variable* var, std::string n = "");
    
    DataSet(std::vector<Variable*>& vars, std::string n = "");
    DataSet(std::set<Variable*>& vars, std::string n = "");
    DataSet(std::initializer_list<Variable*> vars, std::string n = "");
    
    virtual ~DataSet() = default;

    virtual void addEvent() = 0; // Must increment numEventsAdded
    virtual void addWeightedEvent(fptype weight);
    
    
    // Helper to add 1D events
    virtual void addEvent (fptype val) ;

    Variable_v::const_iterator begin() const {
        return variables.begin();
    }
    Variable_v::const_iterator end() const {
        return variables.end();
    }
    
    Variable_v getVariables() const;

    size_t size() const {
        return variables.size();
    }
    
    int numEvents() const {
        return numEventsAdded;
    }

    std::string getName() const {
        return name;
    }

protected:
    std::vector<fptype> getCurrentValues() const;
    size_t indexOfVariable(Variable* var) const;
    size_t numEventsAdded {0};
    
private:
    /// Make a name, does not change the exising name. Called by all constructors.
    void generateName();
    
    std::string name;

protected:
    std::vector<Variable*> variables;
};

