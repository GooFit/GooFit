#ifndef VARIABLE_HH
#define VARIABLE_HH

#include <string>
#include <map>
#include <iostream>
#include <cassert>
#include "goofit/GlobalCudaDefines.h"

struct Indexable {
    Indexable(std::string n, fptype val = 0) : name(n), value(val), index(-1) {}

    int getIndex() const {
        return index;
    }
    std::string name;
    fptype value;
    int index;
};

struct Variable : Indexable {
    // Contains information about a parameter allowed
    // to vary in MINUIT, or an observable passed to a
    // data set. The index can refer either to cudaArray
    // or to an event.

    Variable(std::string n);
    Variable(std::string n, fptype val);
    Variable(std::string n, fptype dn, fptype up);
    Variable(std::string n, fptype v, fptype dn, fptype up);
    Variable(std::string n, fptype v, fptype e, fptype dn, fptype up);
    ~Variable();

    fptype error;
    fptype upperlimit;
    fptype lowerlimit;
    int numbins;
    bool fixed;
    fptype blind;
};

struct Constant : Indexable {
    // This is similar to Variable, but the index points
    // to functorConstants instead of cudaArray.

    Constant(std::string n, fptype val) : Indexable(n, val) {}
    ~Constant() {}
};


#endif
