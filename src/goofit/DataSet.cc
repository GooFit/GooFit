#include "goofit/DataSet.h"
#include "goofit/Error.h"
#include "goofit/Variable.h"

#include <climits>
#include <cstdlib>
#include <utility>

namespace GooFit {

DataSet::DataSet(Observable *var, std::string n)
    : name(std::move(n))
    , variables({var}) {
    generateName();
}

DataSet::DataSet(std::vector<Observable *> &vars, std::string n)
    : name(std::move(n))
    , variables(vars) {
    generateName();
}

DataSet::DataSet(std::set<Observable *> &vars, std::string n)
    : name(std::move(n))
    , variables(std::begin(vars), std::end(vars)) {
    generateName();
}

DataSet::DataSet(std::initializer_list<Observable *> vars, std::string n)
    : name(std::move(n))
    , variables(vars) {
    generateName();
}

void DataSet::addWeightedEvent(fptype) {
    throw GooFit::GeneralError("AddWeightedEvent not implemented for this type of DataSet");
}

std::vector<fptype> DataSet::getCurrentValues() const {
    std::vector<fptype> values;

    for(Observable *v : variables) {
        values.push_back(v->getValue());
    }

    return values;
}

const std::vector<Observable *> &DataSet::getVariables() const { return variables; }

size_t DataSet::indexOfVariable(Observable *var) const {
    for(size_t i = 0; i < variables.size(); ++i)
        if(var == variables[i])
            return i;

    throw GooFit::GeneralError("Invalid variable access into dataset!");
}

void DataSet::generateName() {
    // Create default name as list of variables.
    if(name != "")
        return;

    for(Observable *v : variables) {
        if(v != variables[0])
            name += ", ";
        name += v->getName();
    }
}

void DataSet::checkAllVars() const {
    for(Observable *v : variables) {
        if(!*v)
            throw GooFit::OutOfRange(v->getName(), v->getValue(), v->getLowerLimit(), v->getUpperLimit());
    }
}
} // namespace GooFit
