#include "goofit/DataSet.h"
#include "goofit/Variable.h"
#include "goofit/Error.h"

#include <cstdlib>
#include <climits>

DataSet::DataSet(Variable* var, std::string n)
: name(n), variables({var}) {
    generateName();
}

DataSet::DataSet(std::vector<Variable*>& vars, std::string n)
    : name(n), variables(vars) {
     generateName();
}

DataSet::DataSet(std::set<Variable*>& vars, std::string n)
: name(n), variables(std::begin(vars), std::end(vars)) {
    generateName();
}

DataSet::DataSet(std::initializer_list<Variable*> vars, std::string n) :
    name(n), variables(vars) {
    generateName();
}




void DataSet::addEvent(fptype val) {
    // Helper method to avoid the user having to wrap
    // every one-dimensional event in a vector.
    assert(1 == variables.size());
    variables[0]->value = val;
    addEvent();
}

void DataSet::addWeightedEvent(fptype) {
    throw GooFit::GeneralError("AddWeightedEvent not implemented for this type of DataSet");
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

size_t DataSet::indexOfVariable(Variable* var) const {
    for(size_t i = 0; i < variables.size(); ++i)
        if(var == variables[i])
            return i;

    throw GooFit::GeneralError("Invalid variable access into dataset!");
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
