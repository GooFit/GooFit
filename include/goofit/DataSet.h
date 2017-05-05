#pragma once

#include "goofit/GlobalCudaDefines.h"

#include <vector>
#include <set>
#include <string>
#include <initializer_list>

#include "goofit/Variable.h"
#include "goofit/Error.h"

class DataSet {
public:
    DataSet(Variable* var, std::string n = "");
    
    DataSet(std::vector<Variable*>& vars, std::string n = "");
    DataSet(std::set<Variable*>& vars, std::string n = "");
    DataSet(std::initializer_list<Variable*> vars, std::string n = "");
    
    virtual ~DataSet() = default;

    virtual void addEvent() = 0; // Must increment numEventsAdded
    
    virtual void addWeightedEvent(fptype weight);
    

    /// This is a helper that allows multiple values to be passed in instead of relying on the content of the Variables.
    template<typename... Args>
    void addEvent(fptype value, Args... args) {
        std::vector<fptype> values {value, static_cast<fptype>(args)...};
        
        if(values.size() != variables.size())
            throw GooFit::GeneralError("You must pass the correct number of values ({}) to addEvent", variables.size());
        
        for(size_t i=0; i<values.size(); i++)
            variables[i]->setValue(values[i]);
        addEvent();
    }

    
    const std::vector<Variable*>& getVariables() const;

    size_t numVariables() const {
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

