#pragma once

#include "goofit/GlobalCudaDefines.h"

#include <vector>
#include <set>
#include <string>

class Variable;

typedef std::vector<Variable*>::const_iterator varConstIt;
typedef std::vector<Variable*>::const_reverse_iterator varConstRIt;

class DataSet {
public:
    DataSet(Variable* var, std::string n = "");
    DataSet(std::vector<Variable*>& vars, std::string n = "");
    DataSet(std::set<Variable*>& vars, std::string n = "");
    virtual ~DataSet();

    void addEvent();
    void addEvent(fptype val);
    virtual void addEventVector(std::vector<fptype>& vals, fptype weight = 1) = 0;
    void addWeightedEvent(fptype weight);

    varConstIt varsBegin() const {
        return variables.begin();
    }
    varConstIt varsEnd() const {
        return variables.end();
    }
    void getVariables(std::vector<Variable*>& vars);

    varConstRIt varsRBegin() const {
        return variables.rbegin();
    }
    varConstRIt varsREnd() const {
        return variables.rend();
    }
    int numVariables() const {
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
    int numEventsAdded;

private:
    void generateName();

    std::vector<Variable*> variables;
    std::string name;
};

