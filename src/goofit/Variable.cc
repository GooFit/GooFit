#include "goofit/Variable.h"
#include "goofit/Error.h"

#include <algorithm>

namespace GooFit {

int max_index(const std::vector<Variable *> &vars) {
    const Variable *max_ind_ptr
        = *std::max_element(std::begin(vars), std::end(vars), [](const Variable *a, const Variable *b) {
              return a->getIndex() < b->getIndex();
          });
    return max_ind_ptr->getIndex();
}

int max_fitter_index(const std::vector<Variable *> &vars) {
    const Variable *max_ind_ptr
        = *std::max_element(std::begin(vars), std::end(vars), [](const Variable *a, const Variable *b) {
              return a->getFitterIndex() < b->getFitterIndex();
          });
    return max_ind_ptr->getFitterIndex();
}

std::ostream &operator<<(std::ostream &o, const GooFit::Variable &var) {
    o << var.getName() << ": " << var.getValue() << " +/- " << var.getError();
    if(!var.fixed)
        o << " [" << var.getLowerLimit() << ", " << var.getUpperLimit() << "]";
    if(var.getIndex() >= 0)
        o << " GooFit index: " << var.getIndex();
    if(var.getFitterIndex() >= 0)
        o << " Fitter index: " << var.getFitterIndex();
    if(var.blind != 0)
        o << " Blinded";

    return o;
}

std::istream &operator>>(std::istream &i, GooFit::Variable &var) { return i >> var.value; }

} // namespace GooFit
