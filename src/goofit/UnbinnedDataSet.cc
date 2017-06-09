#include "goofit/UnbinnedDataSet.h"
#include "goofit/Error.h"
#include "goofit/Log.h"
#include "goofit/Variable.h"

namespace GooFit {

// Special constructor for one variable
UnbinnedDataSet::UnbinnedDataSet(Variable *var, std::string n)
    : DataSet(var, n) {
    data.resize(1);
}

UnbinnedDataSet::UnbinnedDataSet(std::vector<Variable *> &vars, std::string n)
    : DataSet(vars, n) {
    data.resize(vars.size());
}

UnbinnedDataSet::UnbinnedDataSet(std::set<Variable *> &vars, std::string n)
    : DataSet(vars, n) {
    data.resize(vars.size());
}

UnbinnedDataSet::UnbinnedDataSet(std::initializer_list<Variable *> vars, std::string n)
    : DataSet(vars, n) {
    data.resize(vars.size());
}

fptype UnbinnedDataSet::getValue(Variable *var, size_t idx) const {
    if(idx >= getNumEvents()) {
        throw GooFit::GeneralError("UnbinnedDataSet: Attepted to find {} in event {} when only {} events exits",
                                   var->getName(),
                                   idx,
                                   getNumEvents());
    }

    size_t var_idx = indexOfVariable(var);

    return data[var_idx].at(idx);
}

void UnbinnedDataSet::loadEvent(size_t idx) {
    size_t i = 0;
    for(Variable *v : variables) {
        v->setValue(data.at(i++).at(idx));
    }
}

void UnbinnedDataSet::setValueForAllEvents(Variable *var) {
    size_t ivar = indexOfVariable(var);
    for(size_t i = 0; i < getNumEvents(); i++) {
        data[ivar][i] = var->getValue();
    }
}

void UnbinnedDataSet::addEvent() {
    checkAllVars();
    size_t i = 0;
    for(Variable *v : variables)
        data.at(i++).push_back(v->getValue());
    numEventsAdded++;
}
} // namespace GooFit
