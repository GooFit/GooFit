#pragma once

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

#include "goofit/GlobalCudaDefines.h"

class Indexable {
public:
    Indexable(std::string n, fptype val = 0) : name(n), value(val) {}
    virtual ~Indexable() {}
    
    /// Get the GooFit index
    int getIndex() const {return index;}
    
    /// Set the GooFit index
    void setIndex(int value) {index = value;}
    
    /// Get the index from the fitter
    int getFitterIndex() const {return fitter_index;}
    
    /// Set the index (should be done by the fitter)
    void setFitterIndex(int value) {fitter_index = value;}

    /// The variable name. Should be unique
    std::string name;
    
    /// The value of the variable
    fptype value;
    
protected:
    /// The goofit index, -1 if unset
    int index {-1};
    
    /// The fitter index, -1 if unset
    int fitter_index {-1};
};

/// Contains information about a parameter allowed
/// to vary in MINUIT, or an observable passed to a
/// data set. The index can refer either to cudaArray
/// or to an event.
class Variable : public Indexable {
public:
    friend std::ostream& operator<< (std::ostream& o, const Variable& var);


    /// This is a constant varaible
    Variable(std::string n, fptype v)
      : Indexable(n, v)
      , error(0.002)
      , upperlimit(v + 0.01)
      , lowerlimit(v - 0.01)
      , fixed(true) {}
    
    /// This is an independent variable
    Variable(std::string n, fptype dn, fptype up)
    : Indexable(n)
    , upperlimit(up)
    , lowerlimit(dn) {}
    
    
    /// This is a normal variable, with value and upper/lower limits
    Variable(std::string n, fptype v, fptype dn, fptype up)
      : Indexable(n, v)
      , error(0.1*(up-dn))
      , upperlimit(up)
    , lowerlimit(dn) {}
    
    /// This is a full varaible with error scale as well
    Variable(std::string n, fptype v, fptype e, fptype dn, fptype up)
      : Indexable(n, v)
      , error(e)
      , upperlimit(up)
      , lowerlimit(dn) {}
    
    
    virtual ~Variable() = default;

    fptype error;
    fptype upperlimit;
    fptype lowerlimit;
    
    int numbins {100};
    
    /// This "fixes" the variable (constant)
    bool fixed {false};
    
    /// A blinding value to add
    fptype blind {0};
    
    /// True if the value was unchanged since the last iteration
    bool unchanged_ {false};
    
    /// Check to see if the value has changed this iteration (always true the first time)
    bool changed() const {return !unchanged_;}
};

/// This is used to track event number for MPI versions.
/// A cast is done to know whether the values need to be fixed.
class CountingVariable : public Variable {
public:

    using Variable::Variable;
    virtual ~CountingVariable() = default;
};

/// This is similar to Variable, but the index points
/// to functorConstants instead of cudaArray.
class Constant : public Indexable {
public:

    Constant(std::string n, fptype val) : Indexable(n, val) {}
    virtual ~Constant() {}
};

inline std::ostream& operator<< (std::ostream& o, const Variable& var) {
    o << var.name << ": " << var.value << " +/- " << var.error;
    if(!var.fixed)
        o << " [" << var.lowerlimit << ", " << var.upperlimit << "]";
    if(var.getIndex() >= 0)
        o << " GooFit index: " << var.getIndex();
    if(var.getFitterIndex() >= 0)
        o << " Fitter index: " << var.getFitterIndex();

    return o;
}


/// Get the max index of a variable from a list
inline int max_index(const std::vector<Variable*> &vars) {
    const Variable* max_ind_ptr = *std::max_element(std::begin(vars),
                                              std::end(vars),
                                              [](const Variable *a, const Variable *b)
                                              {return a->getIndex() < b->getIndex();});
    return max_ind_ptr->getIndex();
}
