#include "goofit/DataSet.h"
#include "goofit/Variable.h"

#include <cstdlib>
#include <climits>

DataSet::DataSet(Variable* var, std::string n)
: name(n), variables({var}) {
    
    if(n == "")
        generateName();
}

DataSet::DataSet(std::vector<Variable*>& vars, std::string n)
    : name(n), variables(vars) {
        
    if(n == "")
        generateName();
}

DataSet::DataSet(std::set<Variable*>& vars, std::string n)
: name(n), variables(std::begin(vars), std::end(vars)) {
    
    if(n == "")
        generateName();
}

DataSet::DataSet(std::initializer_list<Variable*> vars, std::string n) :
    name(n), variables(vars) {

    if(n == "")
        generateName();
}

void DataSet::addEvent() {
    std::vector<fptype> vals = getCurrentValues();
    addEventVector(vals);
}

void DataSet::addWeightedEvent(fptype weight) {
     std::vector<fptype> vals = getCurrentValues();
    addEventVector(vals, weight);
}

void DataSet::addEvent(fptype val) {
    // Helper method to avoid the user having to wrap
    // every one-dimensional event in a vector.
    assert(1 == variables.size());

    std::vector<fptype> helper;
    helper.push_back(val);
    addEventVector(helper);
}

 std::vector<fptype> DataSet::getCurrentValues() const {
     std::vector<fptype> values;

for(Variable* v : variables) {
        values.push_back(v->value);
    }

    return values;
}

Variable_v DataSet::getVariables() const {
    return variables;
}

unsigned int DataSet::indexOfVariable(Variable* var) const {
    for(unsigned int i = 0; i < variables.size(); ++i) {
        if(var != variables[i])
            continue;

        return i;
    }

    std::cout << "Error: Attempted to get index of variable "
              << var->name
              << " in DataSet of ";

    for(Variable* v : variables) {
        std::cout << "\n  " << v->name << std::endl;
    }

    std::cout << "\nAborting." << std::endl;
    assert(false);
    abort();
    return 0;
}

void DataSet::generateName() {
    // Create default name as list of variables.
    if(name != "")
        return;
    
    for(Variable* v : variables) {
        if(v != variables[0])
            name += ", ";
        name += v->name;
    }
}
