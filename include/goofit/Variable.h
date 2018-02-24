#pragma once

#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <goofit/GlobalCudaDefines.h>
#include <goofit/Log.h>

// Declaring friends

namespace GooFit {

class Params;
class Minuit1;

/// This class is specifically for any indexable quanity: Observables or Variables.
/// It contains the idea of a limited quanity (min==max means no limit).
/// It can also take a value.
class Indexable {
  protected:
    /// The variable name. Should be unique
    std::string name;

    /// The value of the variable
    std::shared_ptr<fptype> value;

    /// The lower limit
    std::shared_ptr<fptype> lowerlimit{std::make_shared<fptype>(0)};

    /// The upper limit
    std::shared_ptr<fptype> upperlimit{std::make_shared<fptype>(0)};

    /// The goofit index, -1 if unset
    std::shared_ptr<int> index{std::make_shared<int>(-1)};

  public:
    Indexable(std::string n, fptype val = 0)
        : name(n)
        , value(std::make_shared<fptype>(val)) {}

    Indexable(std::string n, fptype val, fptype lowerlimit, fptype upperlimit)
        : name(n)
        , value(std::make_shared<fptype>(val))
        , lowerlimit(std::make_shared<fptype>(lowerlimit))
        , upperlimit(std::make_shared<fptype>(upperlimit)) {}

    Indexable &operator=(const Indexable &) = delete;

    virtual ~Indexable() {
        if(value.use_count() == 1)
            GOOFIT_DEBUG("Destroying Indexable: {}", name);
    }

    /// Get the GooFit index
    int getIndex() const { return *index; }
    /// Set the GooFit index
    void setIndex(int value) { *index = value; }

    /// Get the name
    const std::string &getName() const { return name; }
    // The name cannot be changed
    // void setName(const std::string &val) { name = val; }

    /// Get the value
    fptype getValue() const { return *value; }
    /// Set the value
    void setValue(fptype val) { *value = val; }

    /// Get the upper limit
    fptype getUpperLimit() const { return *upperlimit; }
    /// Set the upper limit
    void setUpperLimit(fptype val) { *upperlimit = val; }

    /// Get the lower limit
    fptype getLowerLimit() const { return *lowerlimit; }
    /// Set the lower limit
    void setLowerLimit(fptype val) { *lowerlimit = val; }

    /// Support var = 3
    const Indexable &operator=(const fptype &val) {
        setValue(val);
        return *this;
    }

    // Utilities

    /// Support fptype val = var
    operator fptype() const { return getValue(); }

    /// Support for using as key in map - Notice this does NOT sort by index!
    bool operator<(const Indexable &other) const { return value.get() < other.value.get(); }

    /// Support for comparison - only if really the same object, not just the same value
    bool operator==(const Indexable &other) const { return value.get() == other.value.get(); }

    /// Support for comparison - only if really the same object, not just the same value
    bool operator!=(const Indexable &other) const { return value.get() != other.value.get(); }

    /// Check to see if in range
    operator bool() const {
        return getUpperLimit() == getLowerLimit() || (getValue() <= getUpperLimit() && getValue() >= getLowerLimit());
    }
};

/// Special class for observables. Used in DataSets.
class Observable : public Indexable {
    friend std::istream &operator>>(std::istream &o, Observable &var);

  protected:
    /// The number of bins (mostly for BinnedData, or plotting help)
    std::shared_ptr<size_t> numbins{std::make_shared<size_t>(100)};

    /// This identifies this as an EventNumber if copies are made
    bool is_event_number{false};

  public:
    Observable(std::string name, fptype lowerlimit, fptype upperlimit)
        : GooFit::Indexable(name, 0, lowerlimit, upperlimit) {}

    ~Observable() override = default;

    Observable &operator=(const Observable &) = delete;

    /// Set the number of bins
    void setNumBins(size_t num) { *numbins = num; }
    /// Get the number of bins
    size_t getNumBins() const { return *numbins; }

    /// Get the bin size, (upper-lower) / bins
    fptype getBinSize() const { return (getUpperLimit() - getLowerLimit()) / getNumBins(); }

    /// Check to see if this is an event number
    bool isEventNumber() const { return is_event_number; }

    /// Support var = 3
    const Observable &operator=(const fptype &val) {
        setValue(val);
        return *this;
    }
};

/// Contains information about a parameter allowed
/// to vary in MINUIT, or an observable passed to a
/// data set. The index can refer either to cudaArray
/// or to an event.
class Variable : public Indexable {
    friend std::ostream &operator<<(std::ostream &o, const Variable &var);

  public:
    /// This provides a key for some special classes to access blind info (passkey)
    class Key {
        friend Params;
        friend Minuit1;

        /// Private constructor
        Key() = default;
    };

    // These classes can now be duplicated safely. No pointers needed.
    Variable(const Variable &) = default;
    Variable(Variable &&)      = default;

    /// This is a constant variable
    Variable(std::string n, fptype v)
        : Indexable(n, v, v + 0.01, v - 0.01)
        , error(std::make_shared<fptype>(0.002)) {
        *fixed = true;
    }

    /// This is a normal variable, with value and upper/lower limits
    Variable(std::string n, fptype v, fptype dn, fptype up)
        : Indexable(n, v, dn, up)
        , error(std::make_shared<fptype>(0.1 * (up - dn))) {}

    /// This is a full variable with error scale as well (Fix added to allow single interface to fixed/free vars)
    Variable(std::string n, fptype v, fptype e, fptype dn, fptype up, bool fix = false)
        : Indexable(n, v, dn, up)
        , error(std::make_shared<fptype>(e))
        , fixed(std::make_shared<bool>(fix)) {}

    ~Variable() override = default;

    Variable &operator=(const Variable &) = delete;

    /// Get the error
    fptype getError() const { return *error; }
    /// Set the error
    void setError(fptype val) { *error = val; }

    /// Get the index from the fitter
    int getFitterIndex() const { return *fitter_index; }
    /// Set the index (should be done by the fitter)
    void setFitterIndex(int value) { *fitter_index = value; }

    /// Check to see if the value has changed this iteration (always true the first time)
    bool getChanged() const { return *changed_; }

    /// Check to see if this is a constant
    bool IsFixed() const { return *fixed; }

    /// Set the fixedness of a variable
    void setFixed(bool fix) { *fixed = fix; }

    /// Check to see if this has been changed since last iteration
    void setChanged(bool val = true) { *changed_ = val; }

    /// Hides the number; the real value is the result minus this value. Cannot be retreived once set.
    void setBlind(fptype val) { *blind = val; }

    /// Protected by special locked key, only a few classes have access to create a Key
    fptype getBlind(const Key &) { return *blind; }

    /// Support var = 3
    const Variable &operator=(const fptype &val) {
        setValue(val);
        return *this;
    }

  protected:
    /// The error
    std::shared_ptr<fptype> error;

    /// A blinding value to add
    std::shared_ptr<fptype> blind{std::make_shared<fptype>(0)};

    /// True if the value was unchanged since the last iteration
    std::shared_ptr<bool> changed_{std::make_shared<bool>(true)};

    /// This "fixes" the variable (constant)
    std::shared_ptr<bool> fixed{std::make_shared<bool>(false)};

    /// The fitter index, -1 if unset
    std::shared_ptr<int> fitter_index{std::make_shared<int>(-1)};
};

/// This is used to track event number for MPI versions.
/// A cast is done to know whether the values need to be fixed.
/// Ugly hack because this internally stores a floating point number!
class EventNumber : public Observable {
  public:
    static constexpr fptype maxint{1L << std::numeric_limits<fptype>::digits};

    EventNumber(std::string name, fptype min = 0, fptype max = EventNumber::maxint)
        : GooFit::Observable(name, min, max) {
        is_event_number = true;
        *value          = 0; // Start with 0 event number by default
    }

    ~EventNumber() override = default;

    EventNumber &operator=(const EventNumber &) = delete;

    /// Support var = 3
    const EventNumber &operator=(const fptype &val) {
        setValue(val);
        return *this;
    }
};

/// Get the max index of a variable from a list
int max_index(const std::vector<Variable> &vars);
int max_index(const std::vector<Observable> &vars);

/// Get the max fitter index of a variable from a list
int max_fitter_index(const std::vector<Variable> &vars);

/// Nice print of Variable
std::ostream &operator<<(std::ostream &o, const GooFit::Variable &var);

/// Nice print of Observable
std::ostream &operator<<(std::ostream &o, const GooFit::Observable &var);

/// Allow Observable to be read in
std::istream &operator>>(std::istream &i, GooFit::Observable &var);
} // namespace GooFit
