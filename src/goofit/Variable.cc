#include "goofit/Variable.h"
#include <cmath>

Variable::Variable(std::string n)
    : Indexable(n) {
}

Variable::Variable(std::string n, fptype v)
    : Indexable(n, v)
    , error(0.002)
    , upperlimit(v + 0.01)
    , lowerlimit(v - 0.01)
    , fixed(true) {
}

Variable::Variable(std::string n, fptype dn, fptype up)
    : Indexable(n)
    , upperlimit(up)
    , lowerlimit(dn){
}

Variable::Variable(std::string n, fptype v, fptype dn, fptype up)
    : Indexable(n, v)
    , error(0.1*(up-dn))
    , upperlimit(up)
    , lowerlimit(dn) {
}

Variable::Variable(std::string n, fptype v, fptype e, fptype dn, fptype up)
    : Indexable(n, v)
    , error(e)
    , upperlimit(up)
    , lowerlimit(dn) {
}

Variable::~Variable() {
}

CountingVariable::CountingVariable(std::string n) : Variable(n) {
}

CountingVariable::CountingVariable(std::string n, fptype v) : Variable(n, v) {
}

CountingVariable::CountingVariable(std::string n, fptype dn, fptype up) : Variable(n, dn, up) {
}

CountingVariable::CountingVariable(std::string n, fptype v, fptype dn, fptype up) : Variable(n, v, dn, up) {
}

CountingVariable::CountingVariable(std::string n, fptype v, fptype e, fptype dn, fptype up) : Variable(n, v, e, dn,
            up) {
}

CountingVariable::~CountingVariable() {
}


/// Nice print for variables
std::ostream& operator<< (std::ostream& o, const Variable& var) {
    o << var.name << ": " << var.value << " +/- " << var.error;
    if(!var.fixed)
        o << " [" << var.lowerlimit << ", " << var.upperlimit << "]";
    return o;
}
