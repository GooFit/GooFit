#include <goofit/DataSet.h>
#include <goofit/Error.h>
#include <goofit/Log.h>
#include <goofit/Variable.h>

#include <climits>
#include <cstdlib>
#include <utility>

namespace GooFit {

DataSet::DataSet(const Observable &var, std::string n)
    : name(std::move(n))
    , observables({var}) {
    generateName();
}

DataSet::DataSet(const std::vector<Observable> &vars, std::string n)
    : name(std::move(n))
    , observables(vars) {
    generateName();
}

DataSet::DataSet(const std::set<Observable> &vars, std::string n)
    : name(std::move(n))
    , observables(std::begin(vars), std::end(vars)) {
    generateName();
}

DataSet::DataSet(std::initializer_list<Observable> vars, std::string n)
    : name(std::move(n))
    , observables(vars) {
    generateName();
}

void DataSet::addWeightedEvent(fptype) {
    throw GooFit::GeneralError("AddWeightedEvent not implemented for this type of DataSet");
}

auto DataSet::getCurrentValues() const -> std::vector<fptype> {
    std::vector<fptype> values;

    for(const Observable &v : observables) {
        values.push_back(v.getValue());
    }

    return values;
}

auto DataSet::getObservables() const -> const std::vector<Observable> & { return observables; }

auto DataSet::indexOfVariable(const Observable &var) const -> size_t {
    for(size_t i = 0; i < observables.size(); ++i)
        if(var == observables[i])
            return i;

    throw GooFit::GeneralError("Invalid variable access into dataset!");
}

void DataSet::generateName() {
    // Create default name as list of variables.
    if(name != "")
        return;

    for(const Observable &v : observables) {
        if(v != observables[0])
            name += ", ";
        name += v.getName();
    }
}

void DataSet::checkAllVars() const {
    for(const Observable &v : observables) {
        if(!v)
            throw GooFit::OutOfRange(v.getName(), v.getValue(), v.getLowerLimit(), v.getUpperLimit());
    }
}
} // namespace GooFit
