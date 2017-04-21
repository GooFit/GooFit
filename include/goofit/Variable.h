#ifndef VARIABLE_HH
#define VARIABLE_HH

#include <string>
#include <map>
#include <iostream>
#include <cassert>
#include "goofit/GlobalCudaDefines.h"

class Indexable {
public:
    Indexable(std::string n, fptype val = 0) : name(n), value(val), index(-1) {}
    virtual ~Indexable() {}

    inline int getIndex() const {
        return index;
    }

    std::string name;
    fptype value;
    int index;
};

/// Contains information about a parameter allowed
/// to vary in MINUIT, or an observable passed to a
/// data set. The index can refer either to cudaArray
/// or to an event.
class Variable : public Indexable {
public:
    friend std::ostream& operator<< (std::ostream& o, const Variable& var);
    


    Variable(std::string n);
    Variable(std::string n, fptype val);
    Variable(std::string n, fptype dn, fptype up);
    Variable(std::string n, fptype v, fptype dn, fptype up);
    Variable(std::string n, fptype v, fptype e, fptype dn, fptype up);
    virtual ~Variable();

    fptype error;
    fptype upperlimit;
    fptype lowerlimit;
    int numbins {100};
    bool fixed {false};
    fptype blind {0};
};


/// This is used to track event number for MPI versions.
/// A cast is done to know whether the values need to be fixed.
class CountingVariable : public Variable {
public:

    CountingVariable(std::string n);
    CountingVariable(std::string n, fptype val);
    CountingVariable(std::string n, fptype dn, fptype up);
    CountingVariable(std::string n, fptype v, fptype dn, fptype up);
    CountingVariable(std::string n, fptype v, fptype e, fptype dn, fptype up);
    virtual ~CountingVariable();
};

/// This is similar to Variable, but the index points
/// to functorConstants instead of cudaArray.
class Constant : public Indexable {
public:

    Constant(std::string n, fptype val) : Indexable(n, val) {}
    virtual ~Constant() {}
};

std::ostream& operator<< (std::ostream& o, const Variable& var);

#endif
