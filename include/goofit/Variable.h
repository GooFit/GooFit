#pragma once

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <memory>

#include "goofit/GlobalCudaDefines.h"
#include <goofit/Log.h>

// Declaring friends

namespace GooFit {

class Params;
class Minuit1;

class Indexable {
  public:
    Indexable(std::string n, fptype val = 0)
        : name(n)
        , value(std::make_shared<fptype>(val)) {}

    // These classes can now be duplicated safely. No pointers needed.
    Indexable(const Indexable &) = default;
    Indexable &operator=(const Indexable &) = default;
    Indexable(Indexable &&)                 = default;

    virtual ~Indexable() {
        GOOFIT_DEBUG("Destroying Variable: {}", name);
    }

    /// Get the GooFit index
    int getIndex() const { return *index; }
    /// Set the GooFit index
    void setIndex(int value) { *index = value; }

    /// Get the index from the fitter
    int getFitterIndex() const { return *fitter_index; }
    /// Set the index (should be done by the fitter)
    void setFitterIndex(int value) { *fitter_index = value; }

    /// Get the name
    const std::string &getName() const { return name; }

    // The name cannot be changed
    // void setName(const std::string &val) { name = val; }

    /// Get the value
    fptype getValue() const { return *value; }
    /// Set the value
    void setValue(fptype val) { *value = val; }

    // Utilities

    /// Support var = 3
    fptype operator=(const fptype &val) {
        setValue(val);
        return val;
    }

    /// Support fptype val = var
    operator fptype() const { return getValue(); }

    /// Support for using as key in map - Notice this does NOT sort by index!
    bool operator <(const Indexable& other) const {
        return value.get() < other.value.get();
    }
    
    /// Support for comparison
    bool operator ==(const Indexable& other) const {
        return value.get() == other.value.get();
    }
    

  protected:
    /// The variable name. Should be unique
    std::string name;

    /// The value of the variable
    std::shared_ptr<fptype> value;

    /// The goofit index, -1 if unset
    std::shared_ptr<int> index{std::make_shared<int>(-1)};

    /// The fitter index, -1 if unset
    std::shared_ptr<int> fitter_index{std::make_shared<int>(-1)};
};

/// Contains information about a parameter allowed
/// to vary in MINUIT, or an observable passed to a
/// data set. The index can refer either to cudaArray
/// or to an event.
class Variable : public Indexable {
    friend std::ostream &operator<<(std::ostream &o, const Variable &var);
    friend std::istream &operator>>(std::istream &o, Variable &var);

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
    Variable &operator=(const Variable &) = default;
    Variable(Variable &&)                 = default;

    /// Support var = 3
    void operator=(const fptype &val) { setValue(val); }

    /// This is a constant varaible
    Variable(std::string n, fptype v)
        : Indexable(n, v)
        , error(std::make_shared<fptype>(0.002))
        , upperlimit(std::make_shared<fptype>(v + 0.01))
        , lowerlimit(std::make_shared<fptype>(v - 0.01)) {
        *fixed = true;
    }

    /// This is an independent variable
    Variable(std::string n, fptype dn, fptype up)
        : Indexable(n)
        , upperlimit(std::make_shared<fptype>(up))
        , lowerlimit(std::make_shared<fptype>(dn)) {}

    /// This is a normal variable, with value and upper/lower limits
    Variable(std::string n, fptype v, fptype dn, fptype up)
        : Indexable(n, v)
        , error(std::make_shared<fptype>(0.1 * (up - dn)))
        , upperlimit(std::make_shared<fptype>(up))
        , lowerlimit(std::make_shared<fptype>(dn)) {}

    /// This is a full varaible with error scale as well
    Variable(std::string n, fptype v, fptype e, fptype dn, fptype up)
        : Indexable(n, v)
        , error(std::make_shared<fptype>(e))
        , upperlimit(std::make_shared<fptype>(up))
        , lowerlimit(std::make_shared<fptype>(dn)) {}

    ~Variable() override = default;

    /// Get the error
    fptype getError() const { return *error; }
    /// Set the error
    void setError(fptype val) { *error = val; }

    /// Get the upper limit
    fptype getUpperLimit() const { return *upperlimit; }
    /// Set the upper limit
    void setUpperLimit(fptype val) { *upperlimit = val; }

    /// Get the lower limit
    fptype getLowerLimit() const { return *lowerlimit; }
    /// Set the lower limit
    void setLowerLimit(fptype val) { *lowerlimit = val; }

    /// Check to see if the value has changed this iteration (always true the first time)
    bool getChanged() const { return *changed_; }

    /// Set the number of bins
    void setNumBins(size_t num) { *numbins = num; }

    /// Get the number of bins
    size_t getNumBins() const { return *numbins; }

    /// Check to see if this is a constant
    bool IsFixed() const { return *fixed; }

    /// Set the fixedness of a variable
    void setFixed(bool fix) { *fixed = fix; }

    /// Check to see if this has been changed since last iteration
    void setChanged(bool val = true) { *changed_ = val; }

    /// Get the bin size, (upper-lower) / bins
    fptype getBinSize() const { return (getUpperLimit() - getLowerLimit()) / getNumBins(); }

    /// Hides the number; the real value is the result minus this value. Cannot be retreived once set.
    void setBlind(fptype val) { *blind = val; }

    /// Protected by special locked key, only a few classes have access to create a Key
    fptype getBlind(const Key &) { return *blind; }

    /// Check to see if in range
    operator bool() const { return getValue() <= getUpperLimit() && getValue() >= getLowerLimit(); }

  protected:
    /// The error
    std::shared_ptr<fptype> error;

    /// The upper limit
    std::shared_ptr<fptype> upperlimit;

    /// The lower limit
    std::shared_ptr<fptype> lowerlimit;

    /// A blinding value to add
    std::shared_ptr<fptype> blind{std::make_shared<fptype>(0)};

    /// The number of bins (mostly for BinnedData, or plotting help)
    std::shared_ptr<size_t> numbins{std::make_shared<size_t>(100)};

    /// True if the value was unchanged since the last iteration
    std::shared_ptr<bool> changed_{std::make_shared<bool>(true)};

    /// This "fixes" the variable (constant)
    std::shared_ptr<bool> fixed{std::make_shared<bool>(false)};
};

/// This is used to track event number for MPI versions.
/// A cast is done to know whether the values need to be fixed.
/// Ugly hack because this internally stores a floating point number!
class CountingVariable : public Variable {
  public:
    static constexpr fptype maxint{1L << std::numeric_limits<fptype>::digits};

    using Variable::Variable;
    CountingVariable(std::string name)
        : CountingVariable(name, 0, CountingVariable::maxint) {}

    ~CountingVariable() override = default;

    CountingVariable(const CountingVariable &) = default;
    CountingVariable &operator=(const CountingVariable &) = default;
    CountingVariable(CountingVariable &&)                 = default;

    /// Support var = 3
    void operator=(const fptype &val) { setValue(val); }
};

/// This is similar to Variable, but the index points
/// to functorConstants instead of cudaArray.
class Constant : public Indexable {
  public:
    // These classes can not be duplicated
    Constant(const Constant &) = default;
    Constant &operator=(const Constant &) = default;
    Constant(Constant &&)                 = default;

    /// Support var = 3
    void operator=(const fptype &val) { setValue(val); }

    Constant(std::string n, fptype val)
        : Indexable(n, val) {}
    ~Constant() override = default;
};

/// Get the max index of a variable from a list
int max_index(const std::vector<Variable *> &vars);

/// Get the max fitter index of a variable from a list
int max_fitter_index(const std::vector<Variable *> &vars);

/// Nice print of Variable
std::ostream &operator<<(std::ostream &o, const GooFit::Variable &var);

/// Allow Variable to be read in
std::istream &operator>>(std::istream &i, GooFit::Variable &var);
} // namespace GooFit
