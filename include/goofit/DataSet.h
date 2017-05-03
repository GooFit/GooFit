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

    void addEvent();
    void addEvent(fptype val);
    virtual void addEventVector(std::vector<fptype>& vals, fptype weight = 1) = 0;
    void addWeightedEvent(fptype weight);

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
    unsigned int indexOfVariable(Variable* var) const;
    int numEventsAdded {0};
    
private:
    void generateName();
    std::string name;

protected:
    std::vector<Variable*> variables;
};

