#include <goofit/Error.h>
#include <goofit/Log.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <Eigen/Dense>

namespace GooFit {

// Special constructor for one variable
UnbinnedDataSet::UnbinnedDataSet(const Observable &var, std::string n)
    : DataSet(var, n) {
    data.resize(1);
}

UnbinnedDataSet::UnbinnedDataSet(const std::vector<Observable> &vars, std::string n)
    : DataSet(vars, n) {
    data.resize(vars.size());
}

UnbinnedDataSet::UnbinnedDataSet(const std::set<Observable> &vars, std::string n)
    : DataSet(vars, n) {
    data.resize(vars.size());
}

UnbinnedDataSet::UnbinnedDataSet(std::initializer_list<Observable> vars, std::string n)
    : DataSet(vars, n) {
    data.resize(vars.size());
}

auto UnbinnedDataSet::getValue(const Observable &var, size_t idx) const -> fptype {
    if(idx >= getNumEvents()) {
        throw GooFit::GeneralError("UnbinnedDataSet: Attepted to find {} in event {} when only {} events exits",
                                   var.getName(),
                                   idx,
                                   getNumEvents());
    }

    size_t var_idx = indexOfVariable(var);

    return data[var_idx].at(idx);
}

void UnbinnedDataSet::loadEvent(size_t idx) {
    size_t i = 0;
    for(Observable &v : observables) {
        v.setValue(data.at(i++).at(idx));
    }
}

void UnbinnedDataSet::setValueForAllEvents(const Observable &var) {
    size_t ivar = indexOfVariable(var);
    for(size_t i = 0; i < getNumEvents(); i++) {
        data[ivar][i] = var.getValue();
    }
}

void UnbinnedDataSet::addEvent() {
    checkAllVars();
    size_t i = 0;
    for(const Observable &v : observables)
        data.at(i++).push_back(v.getValue());
    numEventsAdded++;
}

// Utility function to make a grid of any dimisinion
void make_a_grid(std::vector<Observable> ret, UnbinnedDataSet *grid, Observable *eventNum) {
    if(ret.empty()) {
        grid->addEvent();

        if(eventNum != nullptr)
            eventNum->setValue(eventNum->getValue() + 1);

        return;
    }

    Observable var = ret.back();
    ret.pop_back(); // safe because this is a copy

    if(var.isEventNumber()) { // Skip
        make_a_grid(ret, grid, eventNum);
    } else {
        for(int i = 0; i < var.getNumBins(); ++i) {
            double step = (var.getUpperLimit() - var.getLowerLimit()) / var.getNumBins();
            var.setValue(var.getLowerLimit() + (i + 0.5) * step);
            make_a_grid(ret, grid, eventNum);
        }
    }
}

void UnbinnedDataSet::fillWithGrid() {
    for(const auto &dat : data)
        if(!dat.empty())
            throw GeneralError("The dataset must be empty to fill with grid!");

    std::vector<Observable> ret = getObservables();

    Observable *evt_num = nullptr;

    for(Observable &val : ret)
        if(val.isEventNumber())
            evt_num = &val;

    make_a_grid(ret, this, evt_num);
}

} // namespace GooFit
