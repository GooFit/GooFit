#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <initializer_list>
#include <set>
#include <string>
#include <vector>

#include <goofit/Error.h>
#include <goofit/Variable.h>

namespace GooFit {

class DataSet {
  public:
    DataSet(const Observable &var, std::string n = "");

    DataSet(const std::vector<Observable> &vars, std::string n = "");
    DataSet(const std::set<Observable> &vars, std::string n = "");
    DataSet(std::initializer_list<Observable> vars, std::string n = "");

    virtual ~DataSet() = default;

    virtual void addEvent() = 0; // Must increment numEventsAdded

    virtual void addWeightedEvent(fptype weight);

    /// This is a helper that allows multiple values to be passed in instead of relying on the content of the Variables.
    template <typename... Args>
    void addEvent(fptype value, Args... args) {
        std::vector<fptype> values{value, static_cast<fptype>(args)...};

        if(values.size() != observables.size())
            throw GooFit::GeneralError("You must pass the correct number of values ({}) to addEvent",
                                       observables.size());

        for(size_t i = 0; i < values.size(); i++)
            observables[i].setValue(values[i]);
        addEvent();
    }

    auto getObservables() const -> const std::vector<Observable> &;

    auto numVariables() const -> size_t { return observables.size(); }

    auto getNumEvents() const -> size_t { return numEventsAdded; }

    auto getName() const -> std::string { return name; }

  protected:
    auto getCurrentValues() const -> std::vector<fptype>;
    auto indexOfVariable(const Observable &var) const -> size_t;
    size_t numEventsAdded{0};

    /// Throw an error if any variables are out of range, call in addEvent
    void checkAllVars() const;

  private:
    /// Make a name, does not change the existing name. Called by all constructors.
    void generateName();

    std::string name;

  protected:
    std::vector<Observable> observables;
};

} // namespace GooFit
